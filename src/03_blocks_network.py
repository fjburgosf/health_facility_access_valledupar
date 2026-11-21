"""Step 03: urban blocks, census population and health questions, and network distances.

- Blocks: urban blocks (class 1) of the DANE IPM layer, with their area-weighted centroid.
- Population: persons of the CNPV 2018 linked to each block through the MGN file (COD_DANE_ANM).
- Persons with a health problem in the last 30 days (P_ENFERMO = 1): formal care = P_QUEHIZO_PPAL in {1, 2}
  (went to their health insurer or to a private doctor); 3-9 = other responses; 11, 99 and the undocumented
  code 0 excluded.
- Road graph: OpenStreetMap ways as an undirected graph with edge lengths in metres. Blocks and facilities are
  snapped to the nearest node of the largest connected component; distances come from a multi-source Dijkstra
  that starts from a virtual node joined to every facility node by its snapping distance.

Outputs: data/processed/blocks.csv, ill_persons.csv; results/audit_blocks.json

Usage (from the repository root):  python src/03_blocks_network.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CNPV = ROOT / "data" / "cnpv2018_cesar_csv"
OUT = ROOT / "data" / "processed"
RES = ROOT / "results"
LAT0, LON0 = 10.46, -73.25
KX, KY = 111320 * np.cos(np.radians(LAT0)), 110540


def xy(lat, lon):
    return np.c_[(np.asarray(lon) - LON0) * KX, (np.asarray(lat) - LAT0) * KY]


def centroid(ring):
    P = xy(ring[:, 1], ring[:, 0])
    x, y = P[:, 0], P[:, 1]
    a = x[:-1] * y[1:] - x[1:] * y[:-1]
    A = a.sum() / 2
    if abs(A) < 1e-9:
        return P.mean(axis=0), 0.0
    return np.array([((x[:-1] + x[1:]) * a).sum() / (6 * A), ((y[:-1] + y[1:]) * a).sum() / (6 * A)]), abs(A)


def blocks():
    fc = json.loads((RAW / "ipm_manzana_valledupar.geojson").read_text(encoding="utf-8"))
    rows = []
    for f in fc["features"]:
        pr = f["properties"]
        if pr["COD_DANE"][5] != "1":
            continue
        g = f["geometry"]
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        cs, ws = [], []
        for poly in polys:
            c, w = centroid(np.array(poly[0]))
            cs.append(c)
            ws.append(w)
        c = np.average(np.array(cs), axis=0, weights=np.array(ws) if sum(ws) > 0 else None)
        rows.append({"block": pr["COD_DANE"], "ipm": pr["ipm"], "x": c[0], "y": c[1], "area_m2": sum(ws)})
    return pd.DataFrame(rows).drop_duplicates("block")


def census():
    mgn = pd.read_csv(CNPV / "CNPV2018_MGN_A2_20.CSV", dtype=str,
                      usecols=["U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA", "COD_DANE_ANM"])
    mgn = mgn[(mgn.U_MPIO == "001") & (mgn.UA_CLASE == "1")]
    per = pd.read_csv(CNPV / "CNPV2018_5PER_A2_20.CSV", dtype=str,
                      usecols=["U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA", "P_SEXO", "P_EDADR", "P_ENFERMO",
                               "P_QUEHIZO_PPAL"])
    per = per[(per.U_MPIO == "001") & (per.UA_CLASE == "1")]
    per = per.merge(mgn[["COD_ENCUESTAS", "U_VIVIENDA", "COD_DANE_ANM"]], on=["COD_ENCUESTAS", "U_VIVIENDA"], how="left")
    # socioeconomic stratum of the dwelling (official utility stratification, unrelated to health; )
    viv = pd.read_csv(CNPV / "CNPV2018_1VIV_A2_20.CSV", dtype=str, usecols=["U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA", "VA1_ESTRATO"])
    viv = viv[(viv.U_MPIO == "001") & (viv.UA_CLASE == "1")]
    per = per.merge(viv[["COD_ENCUESTAS", "U_VIVIENDA", "VA1_ESTRATO"]], on=["COD_ENCUESTAS", "U_VIVIENDA"], how="left")
    return per


def sym_min_csr(I, J, W, n):
    """Symmetric sparse graph keeping the minimum weight of repeated (i, j) pairs; scipy would add them up."""
    I, J, W = np.r_[I, J], np.r_[J, I], np.r_[W, W]
    order = np.lexsort((W, J, I))
    I, J, W = I[order], J[order], W[order]
    first = np.r_[True, (I[1:] != I[:-1]) | (J[1:] != J[:-1])]
    return coo_matrix((W[first], (I[first], J[first])), shape=(n, n)).tocsr()


def graph():
    d = json.loads((RAW / "osm_roads.json").read_text(encoding="utf-8"))
    ids, coords, I, J, W = {}, [], [], [], []
    for e in d["elements"]:
        g = e.get("geometry")
        if not g:
            continue
        P = xy([p["lat"] for p in g], [p["lon"] for p in g])
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
    A = sym_min_csr(np.array(I), np.array(J), np.array(W), n)
    ncomp, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax()
    return A, np.array(coords), lab == big


def network_distance(A, C, keep, sources, targets):
    """Distance from each target point to the nearest source point along the graph (plus snapping distances)."""
    nodes = np.where(keep)[0]
    tree = cKDTree(C[nodes])
    ds, js = tree.query(sources)
    dt, jt = tree.query(targets)
    n = A.shape[0]
    B = A.tocoo()
    # virtual source node n joined to the node of every site; several sites on one node keep the minimum distance
    G = sym_min_csr(np.r_[B.row, np.full(len(js), n)], np.r_[B.col, nodes[js]], np.r_[B.data, ds], n + 1)
    dist = dijkstra(G, directed=False, indices=n)
    return dist[nodes[jt]] + dt, dt


def main():
    OUT.mkdir(exist_ok=True)
    b = blocks()
    per = census()
    pop = per.groupby("COD_DANE_ANM").size().rename("persons")
    b = b.merge(pop, left_on="block", right_index=True, how="left").fillna({"persons": 0})
    s1 = per.assign(s1=(per.VA1_ESTRATO == "1").astype(float)).groupby("COD_DANE_ANM").s1.mean().rename("strat1")
    b = b.merge(s1, left_on="block", right_index=True, how="left")
    fac = pd.read_csv(OUT / "facilities.csv")
    fac = fac[fac.urban]
    A, C, keep = graph()
    P = b[["x", "y"]].values
    for col, F in (("d_ips", fac), ("d_hosp", fac[fac.hospital])):
        S = F[["x", "y"]].values
        b[col], snap = network_distance(A, C, keep, S, P)
        b[col.replace("d_", "e_")] = cKDTree(S).query(P)[0]
    b["snap_m"] = snap
    b.to_csv(OUT / "blocks.csv", index=False)
    ill = per[per.P_ENFERMO == "1"].copy()
    ill = ill[ill.P_QUEHIZO_PPAL.isin([str(k) for k in range(1, 10)])]  # 0 (undocumented), 11 and 99 excluded
    ill["formal"] = ill.P_QUEHIZO_PPAL.isin(["1", "2"]).astype(int)
    ill = ill.rename(columns={"COD_DANE_ANM": "block"})[["block", "P_SEXO", "P_EDADR", "P_QUEHIZO_PPAL", "formal", "VA1_ESTRATO"]]
    ill.to_csv(OUT / "ill_persons.csv", index=False)
    res = b[(b.persons > 0) & b.ipm.notna()]
    audit = {"urban_blocks_ipm_layer": int(len(b)), "residential_blocks": int(len(res)),
             "persons_urban_cnpv": int(len(per)), "persons_linked_to_block": int(per.COD_DANE_ANM.notna().sum()),
             "persons_in_residential_blocks": int(res.persons.sum()),
             "persons_in_blocks_without_ipm_layer": int(len(per) - per.COD_DANE_ANM.isin(b.block).sum()),
             "graph_nodes": int(A.shape[0]), "graph_nodes_largest_component": int(keep.sum()),
             "facilities_urban": int(len(fac)), "hospitals_urban": int(fac.hospital.sum()),
             "snap_m_median": float(res.snap_m.median()), "snap_gt_200": int((res.snap_m > 200).sum()),
             "ill_persons": int((per.P_ENFERMO == "1").sum()), "ill_with_response": int(len(ill)),
             "ill_formal_share": float(ill.formal.mean()),
             "p_quehizo_counts": ill.P_QUEHIZO_PPAL.value_counts().sort_index().to_dict()}
    (RES / "audit_blocks.json").write_text(json.dumps(audit, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=1))
    print(res[["d_ips", "d_hosp", "e_ips", "ipm"]].describe().round(1))


if __name__ == "__main__":
    main()
