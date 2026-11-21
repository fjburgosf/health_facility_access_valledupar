"""Step 01: download the open sources, with provenance.

- IPM by city block, DANE (ArcGIS VULNRB_IPMxMZ, layer 17), Valledupar (20001), polygons in WGS84 (GeoJSON).
- OpenStreetMap through the Overpass API: health facilities and the road network inside the bounding box of the
  urban blocks (only the bounding box is sent).
- REPS (datos.gov.co c36g-9fc2): sites of health providers in Valledupar. Only institutional fields are kept
  (provider class, nature, provider and site names); rows of independent professionals are counted but not stored,
  since their names and addresses are personal data.

The CNPV 2018 microdata of Cesar (ANDA, catalogue 643) are read from data/cnpv2018_cesar_csv/ (CSV files of the persons, dwellings and geographic modules; see README).

Usage (from the repository root):  python src/01_get_sources.py
"""
import datetime as dt
import hashlib
import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
IPM = "https://services5.arcgis.com/PzzAFzyM92dtPsCS/ArcGIS/rest/services/VULNRB_IPMxMZ/FeatureServer/17/query"
OVERPASS = "https://overpass-api.de/api/interpreter"
UA = {"User-Agent": "research-accessibility-study/1.0"}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def get_ipm():
    feats, off = [], 0
    while True:
        r = requests.get(IPM, params={"where": "COD_MPIO='20001'", "outFields": "*", "outSR": 4326, "f": "geojson",
                                      "resultOffset": off, "resultRecordCount": 2000}, headers=UA, timeout=300)
        r.raise_for_status()
        fs = r.json()["features"]
        feats += fs
        if len(fs) < 2000:
            break
        off += 2000
    return {"type": "FeatureCollection", "features": feats}


def bbox(fc):
    xs, ys = [], []
    for f in fc["features"]:
        if f["properties"]["COD_DANE"][5] != "1":
            continue
        g = f["geometry"]
        rings = g["coordinates"] if g["type"] == "Polygon" else [r for p in g["coordinates"] for r in p]
        for ring in rings:
            for x, y in ring:
                xs.append(x)
                ys.append(y)
    pad = 0.01
    return min(ys) - pad, min(xs) - pad, max(ys) + pad, max(xs) + pad


def overpass(q):
    for k in range(5):
        r = requests.post(OVERPASS, data={"data": q}, headers=UA, timeout=600)
        if r.status_code == 200:
            return r.content
        time.sleep(20 * (k + 1))
    r.raise_for_status()


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    prov = {}
    fc = get_ipm()
    b = json.dumps(fc).encode()
    (RAW / "ipm_manzana_valledupar.geojson").write_bytes(b)
    prov["ipm_manzana_valledupar.geojson"] = {"url": IPM, "where": "COD_MPIO='20001'", "features": len(fc["features"]),
                                             "downloaded_utc": now, "sha256": sha(b), "terms": "DANE"}
    s, w, n, e = bbox(fc)
    box = f"({s:.5f},{w:.5f},{n:.5f},{e:.5f})"
    q_health = (f'[out:json][timeout:300];(nwr["amenity"~"^(hospital|clinic|doctors)$"]{box};'
                f'nwr["healthcare"~"^(hospital|clinic|centre|doctor)$"]{box};);out center tags;')
    hb = overpass(q_health)
    (RAW / "osm_health.json").write_bytes(hb)
    prov["osm_health.json"] = {"url": OVERPASS, "query": q_health, "elements": len(json.loads(hb)["elements"]),
                               "downloaded_utc": now, "sha256": sha(hb), "licence": "ODbL (OpenStreetMap contributors)"}
    q_roads = (f'[out:json][timeout:600];way["highway"~"^(motorway|trunk|primary|secondary|tertiary|unclassified|'
               f'residential|living_street|service|pedestrian|track|footway|path|steps|cycleway|road|motorway_link|'
               f'trunk_link|primary_link|secondary_link|tertiary_link)$"]{box};out geom;')
    rb = overpass(q_roads)
    (RAW / "osm_roads.json").write_bytes(rb)
    prov["osm_roads.json"] = {"url": OVERPASS, "query": q_roads, "elements": len(json.loads(rb)["elements"]),
                              "downloaded_utc": now, "sha256": sha(rb), "licence": "ODbL (OpenStreetMap contributors)"}
    r = requests.get("https://www.datos.gov.co/resource/c36g-9fc2.json",
                     params={"$where": "upper(municipiosededesc)='VALLEDUPAR'", "$limit": 50000}, headers=UA, timeout=300)
    r.raise_for_status()
    rows = r.json()
    keep = [{"claseprestador": x.get("claseprestador"), "naturalezajuridica": x.get("naturalezajuridica"),
             "ese": x.get("ese"), "nombreprestador": x.get("nombreprestador"), "nombresede": x.get("nombresede"),
             "codigohabilitacionsede": x.get("codigohabilitacionsede"),
             # institutional address, kept only for IPS (legal entities) and used for local geocoding
             "direccion_sede": x.get("direcci_nsede") if (x.get("claseprestador") or "").startswith("Instituciones") else None,
             "fecha_corte_reps": x.get("fecha_corte_reps")}
            for x in rows if "Profesional" not in (x.get("claseprestador") or "")]
    counts = {}
    for x in rows:
        counts[x.get("claseprestador")] = counts.get(x.get("claseprestador"), 0) + 1
    out = {"counts_by_class": counts, "institutional_sites": keep}
    ob = json.dumps(out, ensure_ascii=False, indent=1).encode("utf-8")
    (RAW / "reps_valledupar_institutional.json").write_bytes(ob)
    prov["reps_valledupar_institutional.json"] = {"url": r.url, "rows_total": len(rows), "rows_kept": len(keep),
                                                  "downloaded_utc": now, "sha256": sha(ob), "licence": "CC BY-SA 4.0"}
    (RAW / "PROVENANCE.json").write_text(json.dumps(prov, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("query",)} for k, v in prov.items()},
                     ensure_ascii=False, indent=1))
    print("counts by class:", counts)


if __name__ == "__main__":
    main()
