[out:csv(::type,::id,::lat,::lon,name,place,"name:es","name:pt";true;"\t")][timeout:900][maxsize:1073741824];
(area["ISO3166-1"="ES"][admin_level=2];area["ISO3166-1"="PT"][admin_level=2];area["ISO3166-1"="AD"][admin_level=2];)->.a;
nwr(area.a)["place"~"^(city|town|village|hamlet|isolated_dwelling|locality|farm)$"]["name"];
out center;
