"""Step 01b: OpenStreetMap road network as it was in mid-2018 (attic query), to match the 2018 census.

Same highway classes and bounding box as step 01, with the Overpass [date:...] setting, which returns the ways and
their geometry as they existed at that date. The output has the same format as osm_roads.json, so the graph of
step 03 can be built from either file.

Output: data/raw/osm_roads_2018.json (and its entry in PROVENANCE.json)

Usage (from the repository root):  python src/01b_get_roads_2018.py
"""
import datetime as dt
import hashlib
import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OVERPASS = "https://overpass-api.de/api/interpreter"
UA = {"User-Agent": "research-accessibility-study/1.0"}
DATE = "2018-07-01T00:00:00Z"   # middle of the census fieldwork year
BOX = "(10.40590,-73.30764,10.51685,-73.20994)"
CLASSES = ("motorway|trunk|primary|secondary|tertiary|unclassified|residential|living_street|service|pedestrian|"
           "track|footway|path|steps|cycleway|road|motorway_link|trunk_link|primary_link|secondary_link|tertiary_link")


def main():
    q = f'[out:json][timeout:900][date:"{DATE}"];way["highway"~"^({CLASSES})$"]{BOX};out geom;'
    for k in range(5):
        r = requests.post(OVERPASS, data={"data": q}, headers=UA, timeout=900)
        if r.status_code == 200:
            break
        time.sleep(30 * (k + 1))
    r.raise_for_status()
    b = r.content
    (RAW / "osm_roads_2018.json").write_bytes(b)
    prov_path = RAW / "PROVENANCE.json"
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    prov["osm_roads_2018.json"] = {"url": OVERPASS, "query": q, "elements": len(json.loads(b)["elements"]),
                                   "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                                   "sha256": hashlib.sha256(b).hexdigest(), "licence": "ODbL (OpenStreetMap contributors)"}
    prov_path.write_text(json.dumps(prov, ensure_ascii=False, indent=1), encoding="utf-8")
    print("ways in 2018:", prov["osm_roads_2018.json"]["elements"])


if __name__ == "__main__":
    main()
