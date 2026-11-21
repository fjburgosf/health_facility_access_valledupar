"""Step 03b: health-facility network and road network contemporary with the 2018 census .

- Opening dates and enabled services of the IPS sites come from the REPS portal (Registro actual, Cesar /
  Valledupar, IPS; read on 2026-11-21): sedes_institucional.csv and servicios_flags_sede.txt.
- 2018 network: urban geocoded sites (step 02) whose opening date is on or before 2018-07-01. Sites that closed
  between 2018 and 2026 are no longer in the current register, so this network is a lower bound of the 2018 supply.
- Primary care: sites with general or family medicine outpatient consultation (REPS services 328 or 325).
- Emergency: sites with an emergency service (1102).
- Road graphs: OpenStreetMap 2026 (step 01) and 2018 (step 01b).
- Geocoding sensitivity: the 2018 sites displaced at random within 100 m and 250 m (five draws each, seed 20261121).

Output: data/processed/blocks_2018.csv (network distances by block for each network and road year),
results/network_2018_audit.json

Usage (from the repository root):  python src/03b_network_2018.py
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PORTAL = RAW / "reps_portal_2026-11-21"
DAT = ROOT / "data" / "processed"
RES = ROOT / "results"
CENSUS_DATE = "20180701"

_spec = importlib.util.spec_from_file_location("a03", Path(__file__).with_name("03_blocks_network.py"))
a03 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a03)


def graph(path):
    d = json.loads(path.read_text(encoding="utf-8"))
    ids, coords, I, J, W = {}, [], [], [], []
    for e in d["elements"]:
        g = e.get("geometry")
        if not g:
            continue
        P = a03.xy([p["lat"] for p in g], [p["lon"] for p in g])
        keys = [(round(p["lat"], 7), round(p["lon"], 7)) for p in g]
        for k, p in zip(keys, P):
            if k not in ids:
                ids[k] = len(coords)
                coords.append(p)
        for (k1, p1), (k2, p2) in zip(zip(keys[:-1], P[:-1]), zip(keys[1:], P[1:])):
            I.append(ids[k1])
            J.append(ids[k2])
            W.append(float(np.hypot(*(p2 - p1))))
    n = len(coords)
    A = coo_matrix((W + W, (I + J, J + I)), shape=(n, n)).tocsr()
    ncomp, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax()
    return A, np.array(coords), lab == big


def services():
    flags = {}
    for line in (PORTAL / "servicios_flags_sede.txt").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        for tok in line.split():
            k, v = tok.split(":")
            code, sede = k.split("-")
            flags["20001" + code + sede] = int(v)
    return flags


def main():
    s = pd.read_csv(PORTAL / "sedes_institucional.csv", dtype=str)
    s["site"] = s.codigo_habilitacion.str.strip() + s.numero_sede.str.strip()
    fl = services()
    f = pd.read_csv(DAT / "facilities.csv", dtype={"site": str})
    f = f.merge(s[["site", "fecha_apertura"]], on="site", how="left")
    f["svc"] = f.site.map(fl)
    f["open2018"] = f.fecha_apertura.notna() & (f.fecha_apertura <= CENSUS_DATE)
    f["gp"] = (f["svc"].fillna(0).astype(int) & 1) > 0
    f["er"] = (f["svc"].fillna(0).astype(int) & 2) > 0
    fu = f[f.urban].copy()
    nets = {"all2026": fu, "all2018": fu[fu.open2018], "gp2026": fu[fu.gp], "gp2018": fu[fu.gp & fu.open2018],
            "er2026": fu[fu.er], "er2018": fu[fu.er & fu.open2018],
            "all2018_undated": fu[fu.open2018 | fu.fecha_apertura.isna()],
            "hosp2026": fu[fu.hospital], "hosp2018": fu[fu.hospital & fu.open2018]}
    # sites of the register that were not geocoded, by zone of the register
    allsites = s.assign(geocoded=s.site.isin(f.site[f.urban]))
    notgeo = allsites[~allsites.geocoded]
    b = pd.read_csv(DAT / "blocks.csv", dtype={"block": str})
    P = b[["x", "y"]].values
    out = b[["block", "persons", "ipm", "d_ips"]].copy()
    audit = {"urban_sites": int(len(fu)), "urban_with_date": int(fu.fecha_apertura.notna().sum()),
             "urban_open_by_2018_07": int(fu.open2018.sum()),
             "urban_opened_after": int((fu.fecha_apertura.notna() & ~fu.open2018).sum()),
             "urban_with_services": int(fu["svc"].notna().sum()), "urban_gp": int(fu.gp.sum()), "urban_er": int(fu.er.sum()),
             "opening_year_counts_urban": fu.fecha_apertura.str[:4].value_counts().sort_index().to_dict(),
             "networks": {k: int(len(v)) for k, v in nets.items()}, "roads": {},
             "register_sites": int(len(allsites)), "register_not_geocoded": int(len(notgeo)),
             "register_not_geocoded_by_zone": notgeo.tipo_zona.str.strip().value_counts().to_dict()}
    rng = np.random.default_rng(20261121)
    F18 = nets["all2018"][["x", "y"]].values
    jitter = {}
    for radius in (100, 250):
        for k in range(5):
            ang = rng.uniform(0, 2 * np.pi, len(F18))
            rad = radius * np.sqrt(rng.uniform(0, 1, len(F18)))
            jitter[f"jit{radius}_{k}"] = F18 + np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
    for roads, path in (("osm2026", RAW / "osm_roads.json"), ("osm2018", RAW / "osm_roads_2018.json")):
        A, C, keep = graph(path)
        audit["roads"][roads] = {"nodes": int(A.shape[0]), "nodes_largest_component": int(keep.sum())}
        for name, F in nets.items():
            d, snap = a03.network_distance(A, C, keep, F[["x", "y"]].values, P)
            out[f"d_{name}_{roads}"] = d
            if name == "all2018":
                out[f"snap_{roads}"] = snap
        out[f"e_all2018"] = cKDTree(F18).query(P)[0]
        node_xy = C[keep]
        audit["roads"][roads]["facility_snap_median_m"] = float(np.median(cKDTree(node_xy).query(F18)[0]))
        audit["roads"][roads]["facility_snap_p90_m"] = float(np.quantile(cKDTree(node_xy).query(F18)[0], 0.9))
        if roads == "osm2018":
            for name, J in jitter.items():
                out[f"d_{name}_{roads}"] = a03.network_distance(A, C, keep, J, P)[0]
    out.to_csv(DAT / "blocks_2018.csv", index=False)
    res = out[(out.persons > 0) & out.ipm.notna()]
    w = res.persons
    audit["mean_distance_weighted"] = {c: float(np.average(res[c], weights=w)) for c in res.columns if c.startswith("d_")}
    audit["snap_median_blocks"] = {r: float(res[f"snap_{r}"].median()) for r in ("osm2026", "osm2018")}
    (RES / "network_2018_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
