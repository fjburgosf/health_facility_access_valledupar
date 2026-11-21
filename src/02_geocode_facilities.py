"""Step 02: local geocoding of the IPS sites of Valledupar.

A Colombian grid address such as "CALLE 16 # 14-90" lies on Calle 16 next to its intersection with Carrera 14.
Each address is parsed into a main street (type, number, letter) and a cross street (number, letter); the site is
placed at the closest point between the OpenStreetMap ways that carry those two names (accepted when the ways come
within 60 m of each other). No address is sent to any external service. Expected positional error: about one block.

Output: data/processed/facilities.csv, results/geocoding_audit.json

Usage (from the repository root):  python src/02_geocode_facilities.py
"""
import json
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
RES = ROOT / "results"
LAT0, LON0 = 10.46, -73.25
KX, KY = 111320 * np.cos(np.radians(LAT0)), 110540


def xy(lat, lon):
    return np.c_[(np.asarray(lon) - LON0) * KX, (np.asarray(lat) - LAT0) * KY]


def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().upper()
    return " ".join(s.replace("º", " ").replace("°", " ").split())


TYPES = {"CALLE": "CALLE", "CALL": "CALLE", "CLL": "CALLE", "CL": "CALLE", "CARRERA": "CARRERA", "CRA": "CARRERA",
         "CR": "CARRERA", "KRA": "CARRERA", "KR": "CARRERA", "CRR": "CARRERA", "CARR": "CARRERA", "K": "CARRERA",
         "TRANSVERSAL": "TRANSVERSAL", "TRANSV": "TRANSVERSAL", "TRV": "TRANSVERSAL", "TV": "TRANSVERSAL",
         "DIAGONAL": "DIAGONAL", "DG": "DIAGONAL", "DIAG": "DIAGONAL"}
CROSS = {"CALLE": ("CARRERA", "TRANSVERSAL"), "DIAGONAL": ("CARRERA", "TRANSVERSAL"),
         "CARRERA": ("CALLE", "DIAGONAL"), "TRANSVERSAL": ("CALLE", "DIAGONAL")}
ADDR = re.compile(r"^(?P<t>[A-Z]+)\.?\s*-?\s*(?P<n>\d{1,3})\s*(?P<l>[A-H])?(?:\s*BIS)?\b.*?"
                  r"(?:#|\bNO\b|\bN\b|\bNUMERO\b|\bNA\b|\bNRO\b|\s)\s*\.?\s*(?P<cn>\d{1,3})\s*(?P<cl>[A-H])?\b")


def parse(addr):
    a = norm(addr)
    a = re.sub(r"N[O°ºª]\.?(?=\s*\d)", " # ", a)
    m = ADDR.match(a)
    if not m or m.group("t") not in TYPES:
        return None
    return TYPES[m.group("t")], int(m.group("n")), m.group("l") or "", int(m.group("cn")), m.group("cl") or ""


def street_index():
    d = json.loads((RAW / "osm_roads.json").read_text(encoding="utf-8"))
    idx = {}
    pat = re.compile(r"^(CALLE|CARRERA|TRANSVERSAL|DIAGONAL)\s+(\d{1,3})\s*([A-H])?\b")
    for e in d["elements"]:
        m = pat.match(norm(e.get("tags", {}).get("name", "")))
        if not m or "geometry" not in e:
            continue
        pts = xy([g["lat"] for g in e["geometry"]], [g["lon"] for g in e["geometry"]])
        for key in ((m.group(1), int(m.group(2)), m.group(3) or ""), (m.group(1), int(m.group(2)), "*")):
            idx.setdefault(key, []).append(pts)
    return {k: np.vstack(v) for k, v in idx.items()}


def densify(P, step=10.0):
    """Points every `step` metres along consecutive vertices (ways are stacked, so long jumps are skipped)."""
    out = [P[:1]]
    for a, b in zip(P[:-1], P[1:]):
        d = np.hypot(*(b - a))
        if 0 < d < 400:
            k = int(d // step)
            if k:
                out.append(a + (b - a) * (np.arange(1, k + 1)[:, None] / (k + 1)))
        out.append(b[None, :])
    return np.vstack(out)


def locate(p, S):
    t, n, l, cn, cl = p
    mains = [(t, n, l)] + ([(t, n, "*")] if l else []) + ([(t, n, "*")] if not l else [])
    for mk in mains:
        if mk not in S:
            continue
        for ct in CROSS[t]:
            for ck, quality in (((ct, cn, cl), "exact"), ((ct, cn, "*"), "number only")):
                if ck not in S:
                    continue
                A, B = densify(S[mk]), densify(S[ck])
                dist, j = cKDTree(B).query(A)
                i = int(np.argmin(dist))
                if dist[i] <= 60:
                    q = quality if mk[2] != "*" or not l else "number only"
                    return (A[i] + B[j[i]]) / 2, q, float(dist[i])
    return None


def main():
    OUT.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    reps = json.loads((RAW / "reps_valledupar_institutional.json").read_text(encoding="utf-8"))
    sites = pd.DataFrame([x for x in reps["institutional_sites"] if (x["claseprestador"] or "").startswith("Instituciones")])
    S = street_index()
    rows, status = [], {}
    for r in sites.itertuples():
        p = parse(r.direccion_sede or "")
        res = locate(p, S) if p else None
        st = "parsed, not found" if p and res is None else ("not parsed" if p is None else res[1])
        status[st] = status.get(st, 0) + 1
        x = y = np.nan
        if res is not None:
            x, y = res[0]
        name = norm(f"{r.nombreprestador} {r.nombresede}")
        rows.append({"site": r.codigohabilitacionsede, "provider": r.nombreprestador, "site_name": r.nombresede,
                     "nature": r.naturalezajuridica, "ese": r.ese, "x": x, "y": y, "geocode": st,
                     "hospital": bool(re.search(r"HOSPITAL|CLINICA|CLINCA", name))})
    f = pd.DataFrame(rows)
    # urban sites: within 300 m of an urban block of the IPM layer
    fc = json.loads((RAW / "ipm_manzana_valledupar.geojson").read_text(encoding="utf-8"))
    cents = []
    for ft in fc["features"]:
        if ft["properties"]["COD_DANE"][5] != "1":
            continue
        g = ft["geometry"]
        ring = np.array(g["coordinates"][0] if g["type"] == "Polygon" else g["coordinates"][0][0])
        cents.append(xy(ring[:, 1].mean(), ring[:, 0].mean())[0])
    tree = cKDTree(np.array(cents))
    ok = f.x.notna()
    f["urban"] = False
    f.loc[ok, "urban"] = tree.query(f.loc[ok, ["x", "y"]].values)[0] <= 300
    f.to_csv(OUT / "facilities.csv", index=False)
    # check against OSM health objects with a name that matches a geocoded site
    osm = json.loads((RAW / "osm_health.json").read_text(encoding="utf-8"))["elements"]
    checks = []
    for e in osm:
        nm = norm(e.get("tags", {}).get("name", ""))
        lat = e.get("lat", e.get("center", {}).get("lat"))
        lon = e.get("lon", e.get("center", {}).get("lon"))
        if len(nm) < 8 or lat is None:
            continue
        key = [w for w in nm.split() if len(w) > 4 and w not in ("HOSPITAL", "CLINICA", "CENTRO", "SALUD")]
        if not key:
            continue
        cand = f[f.urban & (f.provider.map(norm).str.contains(key[0]) | f.site_name.map(norm).str.contains(key[0]))]
        if len(cand):
            p = xy(lat, lon)[0]
            d = np.hypot(cand.x - p[0], cand.y - p[1]).min()
            checks.append({"osm_name": e["tags"]["name"], "keyword": key[0], "nearest_geocoded_m": round(float(d))})
    audit = {"ips_sites": int(len(f)), "status": status, "geocoded": int(ok.sum()), "urban_geocoded": int(f.urban.sum()),
             "urban_hospitals": int((f.urban & f.hospital).sum()), "distinct_urban_points": int(f[f.urban][["x", "y"]].round(0).drop_duplicates().shape[0]),
             "osm_checks": checks,
             "osm_check_median_m": float(np.median([c["nearest_geocoded_m"] for c in checks])) if checks else None}
    (RES / "geocoding_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
