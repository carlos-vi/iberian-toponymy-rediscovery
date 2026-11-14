"""Runs the whole analysis, from the raw data to the figures and the supplementary material.

    python src/run_all.py              # everything except the downloads (uses the snapshot in data/raw/)
    python src/run_all.py --download   # downloads the raw data first (overwrites nothing that already exists
                                       # unless --refresh is also given)
    python src/run_all.py --from 30    # starts at the step whose script number is ≥ 30

Each step is a separate process; its output goes to logs/<script>.log. Total time on a laptop: 35–60 minutes.
"""
import os, sys, subprocess, time
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SRC = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SRC)
LOGS = os.path.join(ROOT, 'logs')
KGRID = '6,8,10,12,14,16,18,20,22,24'

DOWNLOAD = [('01_download_osm.py', []), ('02_download_geonames.py', []), ('03_download_ign.py', []),
            ('04_download_lexicons.py', []), ('05_download_naturalearth.py', [])]
STEPS = [('05_download_naturalearth.py', []),   # base-map layers, fetched once if they are not there
         ('10_build_corpora.py', []), ('11_basemap.py', []),
         ('21_layers.py', ['S', KGRID]), ('21_layers.py', ['SI', KGRID]), ('21_layers.py', ['M', '18']),
         ('21_layers.py', ['I', '18']), ('21_layers.py', ['G', '18']), ('22_layers_robustness.py', []),
         ('30_scan.py', ['SIo', 'SI', 'Go', 'Mo', 'Io']), ('31_scan_analysis.py', []),
         ('40_wedge.py', []), ('41_ena.py', []), ('42_arabic_frontier.py', []), ('43_aque_eque.py', []),
         ('50_figures.py', []), ('60_numbers.py', []), ('70_supplement.py', []),
         (os.path.join('..', 'supplement', 'build.py'), [])]


def run(script, args, extra=()):
    os.makedirs(LOGS, exist_ok=True)
    name = os.path.splitext(os.path.basename(script))[0] + ('_' + args[0] if args and script.startswith('21') else '')
    log = os.path.join(LOGS, name + '.log')
    t0 = time.time()
    print(f'-> {script} {" ".join(args)}', end='  ', flush=True)
    with open(log, 'w', encoding='utf-8') as f:
        r = subprocess.run([sys.executable, '-u', os.path.join(SRC, script), *args, *extra], cwd=SRC, stdout=f,
                           stderr=subprocess.STDOUT, env={**os.environ, 'PYTHONIOENCODING': 'utf-8'})
    print(f'{"ok" if r.returncode == 0 else "FAILED"} ({time.time() - t0:.0f} s; log: {os.path.relpath(log, ROOT)})', flush=True)
    if r.returncode:
        sys.exit(r.returncode)


def main():
    argv = sys.argv[1:]
    start = int(argv[argv.index('--from') + 1]) if '--from' in argv else 0
    if '--download' in argv:
        for s, a in DOWNLOAD:
            run(s, a, ['--refresh'] if '--refresh' in argv else [])
    for s, a in STEPS:
        num = os.path.basename(s)[:2]
        if num.isdigit() and int(num) < start:
            continue
        run(s, a)


if __name__ == '__main__':
    main()
