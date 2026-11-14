"""Converts supplementary.md into a self-contained HTML file (figures embedded) and, when Chrome, Chromium or Edge is
available, into PDF with the browser's headless mode.

The Markdown understood is deliberately simple: headings, paragraphs, bulleted and numbered lists, tables, images
(`![alt|width](path)`), emphasis, inline code, links, sub- and superscripts (~x~, ^x^).
Usage: python supplement/build.py
"""
import os, re, sys, html, shutil, subprocess, base64

HERE = os.path.dirname(os.path.abspath(__file__))
CSS = """
@page { size: A4; margin: 2.1cm 2.2cm; }
body { font-family: Cambria, Georgia, 'Times New Roman', serif; font-size: 10.6pt; line-height: 1.45; color: #1d1d1b;
       max-width: 46em; margin: 2em auto; padding: 0 1em; background: #fff; }
h1 { font-size: 17pt; line-height: 1.25; margin: 0 0 .5em; }
h2 { font-size: 13pt; margin: 1.5em 0 .4em; border-bottom: 1px solid #ddd; padding-bottom: .15em; }
h3 { font-size: 11.2pt; margin: 1.2em 0 .3em; }
p { margin: .45em 0; text-align: justify; hyphens: auto; }
table { border-collapse: collapse; margin: .7em 0 1em; font-family: Calibri, 'Segoe UI', sans-serif; font-size: 8.6pt;
        width: 100%; page-break-inside: avoid; }
th, td { border-bottom: 1px solid #ddd; padding: 3px 5px; text-align: left; vertical-align: top; }
th { border-bottom: 1.5px solid #888; font-weight: 600; }
img { max-width: 100%; display: block; margin: .8em auto .2em; page-break-inside: avoid; }
code { font-family: Consolas, 'DejaVu Sans Mono', monospace; font-size: 9pt; }
hr { border: 0; border-top: 1px solid #ccc; margin: 1.2em 0; }
li { margin: .2em 0; }
ul, ol { padding-left: 1.4em; }
.keep { page-break-inside: avoid; break-inside: avoid; }
"""


def inline(t):
    t = html.escape(t, quote=False)
    # keep numbers with their units, and 'S + I', on one line
    t = re.sub(r'(\d) (%|‰|km\b)', '\\1\u00a0\\2', t)
    t = t.replace('S + I', 'S\u00a0+\u00a0I')
    t = re.sub(r'`([^`]+)`', r'<code>\1</code>', t)
    t = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', t)
    t = re.sub(r'(?<![A-Za-z0-9*])\*(?!\s)(.+?)(?<!\s)\*(?![A-Za-z0-9*])', r'<i>\1</i>', t)
    t = re.sub(r'<i>(-[^<\s]*)</i>', r'<i style="white-space:nowrap">\1</i>', t)   # no line break inside '-ena'
    t = re.sub(r'~([^~\s]+)~', r'<sub>\1</sub>', t)
    t = re.sub(r'\^([^\^\s]+)\^', r'<sup>\1</sup>', t)
    t = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)', r'<a href="\2">\1</a>', t)
    return t


def img_tag(alt, src):
    """![Text|60%](path): the optional suffix |NN% sets the width of the figure."""
    width = '100%'
    if '|' in alt:
        alt, width = alt.split('|', 1)
    path = os.path.normpath(os.path.join(HERE, src))
    data = base64.b64encode(open(path, 'rb').read()).decode()
    return f'<img alt="{html.escape(alt)}" style="width:{width}" src="data:image/png;base64,{data}">'


def convert(md):
    out, i, lines = [], 0, md.split('\n')
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1; continue
        m = re.match(r'^(#{1,3})\s+(.*)$', ln)
        if m:
            n = len(m.group(1)); out.append(f'<h{n}>{inline(m.group(2))}</h{n}>'); i += 1; continue
        if ln.strip() == '---':
            out.append('<hr>'); i += 1; continue
        m = re.match(r'^!\[([^\]]*)\]\(([^)]+)\)\s*$', ln)
        if m:
            out.append(img_tag(m.group(1), m.group(2))); i += 1; continue
        if ln.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].startswith('|'):
                rows.append([c.strip() for c in lines[i].strip().strip('|').split('|')]); i += 1
            head, body = rows[0], rows[2:]
            t = '<table><thead><tr>' + ''.join(f'<th>{inline(c)}</th>' for c in head) + '</tr></thead><tbody>'
            t += ''.join('<tr>' + ''.join(f'<td>{inline(c)}</td>' for c in r) + '</tr>' for r in body) + '</tbody></table>'
            out.append(t); continue
        if re.match(r'^\s*[-*]\s+', ln) or re.match(r'^\s*\d+\.\s+', ln):
            ordered = bool(re.match(r'^\s*\d+\.\s+', ln)); items = []
            while i < len(lines) and (re.match(r'^\s*[-*]\s+', lines[i]) or re.match(r'^\s*\d+\.\s+', lines[i])):
                items.append(re.sub(r'^\s*([-*]|\d+\.)\s+', '', lines[i])); i += 1
            tag = 'ol' if ordered else 'ul'
            out.append(f'<{tag}>' + ''.join(f'<li>{inline(x)}</li>' for x in items) + f'</{tag}>'); continue
        para = [ln]; i += 1
        while i < len(lines) and lines[i].strip() and not re.match(r'^(#{1,3}\s|\||!\[|---$|\s*[-*]\s|\s*\d+\.\s)', lines[i]):
            para.append(lines[i]); i += 1
        out.append(f'<p>{inline(" ".join(para))}</p>')
    # keep each table with its caption, and each figure with its caption, on the same page
    merged, i = [], 0
    while i < len(out):
        a, b = out[i], out[i + 1] if i + 1 < len(out) else ''
        if (a.startswith('<p><b>Table') and b.startswith('<table')) or (a.startswith('<img') and b.startswith('<p><b>Figure')):
            merged.append(f'<div class="keep">{a}{b}</div>'); i += 2
        else:
            merged.append(a); i += 1
    return '\n'.join(merged)


def browser():
    cands = [r'C:\Program Files\Google\Chrome\Application\chrome.exe',
             r'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe',
             r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
             '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome']
    cands += [shutil.which(n) for n in ('google-chrome', 'chromium', 'chromium-browser', 'msedge')]
    return next((c for c in cands if c and os.path.exists(c)), None)


def build(name):
    md = open(os.path.join(HERE, f'{name}.md'), encoding='utf-8').read()
    title = re.search(r'^#\s+(.*)$', md, re.M).group(1)
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>'
            f'<style>{CSS}</style></head><body>{convert(md)}</body></html>')
    hp = os.path.join(HERE, f'{name}.html'); open(hp, 'w', encoding='utf-8').write(page)
    words = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'’.,–-]*", re.sub(r'!\[[^\]]*\]\([^)]*\)', '', md)))
    print(f'{name}: HTML {hp} ({words:,} words in the Markdown source)')
    exe = browser()
    if exe:
        pdf = os.path.join(HERE, f'{name}.pdf')
        url = 'file:///' + hp.replace(os.sep, '/').lstrip('/')
        subprocess.run([exe, '--headless=new', '--disable-gpu', '--no-pdf-header-footer', f'--print-to-pdf={pdf}', url],
                       check=False, capture_output=True, timeout=300)
        print(f'{name}: PDF {pdf} {"written" if os.path.exists(pdf) else "NOT written"}')
    else:
        print('no Chrome/Chromium/Edge found: PDF skipped (print the HTML from a browser instead)')


if __name__ == '__main__':
    for n in (sys.argv[1:] or ['supplementary']):
        if os.path.exists(os.path.join(HERE, f'{n}.md')):
            build(n)
