"""Step 04: co-accumulation of block poverty and distance, and formal care-seeking.

1. Population-weighted distance by IPM quintile; Spearman correlation between IPM and distance.
2. Global Moran's I of IPM and distance, and bivariate local Moran (IPM of the block against the spatial lag of the
   distance), k nearest neighbours (row-standardised), 999 conditional permutations, p < 0.05.
3. Persons with a health problem: logistic regression of formal care on network distance (km), block IPM (per 10
   points), sex and age group, with standard errors clustered by block; crude shares by tertile of distance.
4. Sensitivity: k = 4 and 12; distance to the nearest hospital; Euclidean distance; blocks snapped at > 200 m excluded.

Outputs (results): quintiles.csv, lisa.csv, care_models.csv, care_tertiles.csv, summary.json

Usage (from the repository root):  python src/04_analysis.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DAT = ROOT / "data" / "processed"
RES = ROOT / "results"
SEED = 20261120


def knn_weights(P, k):
    _, idx = cKDTree(P).query(P, k=k + 1)
    return idx[:, 1:]


def z(v):
    return (v - v.mean()) / v.std()


def moran(v, nb):
    zv = z(v)
    return float((zv * zv[nb].mean(axis=1)).sum() / (zv * zv).sum())


def bivariate_lisa(x, y, nb, perms=999, seed=SEED):
    """Local Moran of x against the spatial lag of y; conditional permutation pseudo p-values (two-sided), with the
    k neighbours of each permutation drawn without replacement."""
    rng = np.random.default_rng(seed)
    zx, zy = z(x), z(y)
    k = nb.shape[1]
    I = zx * zy[nb].mean(axis=1)
    n = len(x)
    count = np.zeros(n)
    for i in range(n):
        others = np.delete(np.arange(n), i)
        draws = rng.choice(others, size=(perms, k), replace=True)
        dup = np.array([len(set(row)) < k for row in draws])
        while dup.any():  # each permutation draws k distinct blocks (sampling without replacement)
            draws[dup] = rng.choice(others, size=(int(dup.sum()), k), replace=True)
            dup = np.array([len(set(row)) < k for row in draws])
        Ip = zx[i] * zy[draws].mean(axis=1)
        count[i] = (np.abs(Ip) >= abs(I[i])).sum()
    p = (count + 1) / (perms + 1)
    lag = zy[nb].mean(axis=1)
    q = np.where(zx > 0, np.where(lag > 0, "High-High", "High-Low"), np.where(lag > 0, "Low-High", "Low-Low"))
    return I, p, q


def perm_p_global(stat, x, y, nb, perms=999, seed=SEED):
    """Pseudo p-value of a global statistic under random permutation of y (two-sided)."""
    rng = np.random.default_rng(seed)
    obs = stat(x, y, nb)
    sims = np.array([stat(x, rng.permutation(y), nb) for _ in range(perms)])
    return float(((np.abs(sims) >= abs(obs)).sum() + 1) / (perms + 1))


def bh(p, q=0.05):
    """Benjamini-Hochberg: boolean array of discoveries at false discovery rate q."""
    p = np.asarray(p)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    thr = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.zeros(len(p), bool)
    out[order] = thr <= q
    return out


def load_blocks():
    """Residential blocks with the main distances contemporary with the 2018 census (step 03b):
    IPS sites open by mid-2018 on the 2018 OpenStreetMap roads. The 2026 distances are kept with suffix 2026."""
    b = pd.read_csv(DAT / "blocks.csv", dtype={"block": str})
    b18 = pd.read_csv(DAT / "blocks_2018.csv", dtype={"block": str}).drop(columns=["persons", "ipm", "d_ips"])
    b = b.merge(b18, on="block", how="left")
    b = b.rename(columns={"d_ips": "d_ips2026", "d_hosp": "d_hosp2026", "e_ips": "e_ips2026", "snap_m": "snap_m2026"})
    b["d_ips"], b["d_hosp"] = b.d_all2018_osm2018, b.d_hosp2018_osm2018
    b["e_ips"], b["snap_m"] = b.e_all2018, b.snap_osm2018
    return b[(b.persons > 0) & b.ipm.notna()].reset_index(drop=True)


def global_bivariate(x, y, nb):
    zx, zy = z(x), z(y)
    return float((zx * zy[nb].mean(axis=1)).mean())


def care_model(ill, b, dcol, adjust="ipm", age5=False, cluster_cut=0):
    d = ill.merge(b[["block", dcol, "ipm"]], on="block", how="inner")
    d["dist_km"] = d[dcol] / 1000
    d["ipm10"] = d.ipm / 10
    d["female"] = (d.P_SEXO == "2").astype(int)
    age = pd.to_numeric(d.P_EDADR)
    d["age_grp"] = pd.cut(age, [0, 3, 9, 13, 21], labels=["0-14", "15-44", "45-64", "65+"])
    if age5:  # five-year age groups (P_EDADR codes), 80 or more grouped
        d["age_grp"] = pd.Categorical(np.minimum(age, 17).astype(int))
    est = d.VA1_ESTRATO.fillna("9").map({"1": "1", "2": "2", "3": "3", "4": "4-6", "5": "4-6", "6": "4-6"}).fillna("unknown")
    d["stratum"] = pd.Categorical(est, categories=["1", "2", "3", "4-6", "unknown"])
    cols = ["dist_km", "ipm10", "female", "age_grp"] if adjust == "ipm" else ["dist_km", "female", "age_grp", "stratum"]
    X = pd.get_dummies(d[cols], columns=[c for c in ("age_grp", "stratum") if c in cols], drop_first=True, dtype=float)
    X = sm.add_constant(X)
    groups = pd.factorize(d.block.str[:len(d.block.iloc[0]) - cluster_cut] if cluster_cut else d.block)[0]
    m = sm.GLM(d.formal, X, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": groups})
    ci = m.conf_int()
    out = pd.DataFrame({"or": np.exp(m.params), "lo": np.exp(ci[0]), "hi": np.exp(ci[1]), "p": m.pvalues})
    out.attrs["model"], out.attrs["X"] = m, X
    return out, len(d), d


def main():
    b = load_blocks()
    ill = pd.read_csv(DAT / "ill_persons.csv", dtype={"block": str, "P_SEXO": str, "P_EDADR": str, "VA1_ESTRATO": str})
    P = b[["x", "y"]].values
    summ = {"n_blocks": int(len(b)), "persons": int(b.persons.sum())}
    # 1. quintiles
    b["ipm_q"] = pd.qcut(b.ipm.rank(method="first"), 5, labels=[1, 2, 3, 4, 5])
    wavg = lambda g, c: np.average(g[c], weights=g.persons)
    q = b.groupby("ipm_q", observed=True).apply(lambda g: pd.Series({
        "blocks": len(g), "persons": g.persons.sum(), "ipm_min": g.ipm.min(), "ipm_max": g.ipm.max(),
        "ipm_mean_w": wavg(g, "ipm"), "d_ips_w": wavg(g, "d_ips"), "d_hosp_w": wavg(g, "d_hosp"),
        "share_gt1km": np.average(g.d_ips > 1000, weights=g.persons),
        "d_ips2026_w": wavg(g, "d_ips2026"),
        "d_s18_r26_w": wavg(g, "d_all2018_osm2026"), "d_s26_r18_w": wavg(g, "d_all2026_osm2018"),
        "share_hosp_gt2km": np.average(g.d_hosp > 2000, weights=g.persons)}), include_groups=False)
    q.to_csv(RES / "quintiles.csv")
    for col in ("d_ips", "d_hosp", "e_ips", "d_ips2026", "d_gp2018_osm2018", "d_all2018_osm2026"):
        rho, p = spearmanr(b.ipm, b[col])
        summ[f"spearman_{col}"] = {"rho": float(rho), "p": float(p)}
    # sensitivity: block share of persons in stratum 1, a deprivation measure without health items
    rho, p = spearmanr(b.strat1, b.d_ips)
    summ["spearman_strat1_d_ips"] = {"rho": float(rho), "p": float(p)}
    summ["bimoran_strat1_k8"] = global_bivariate(b.strat1.values, b.d_ips.values, knn_weights(P, 8))
    summ["d_ips_w_all"] = float(wavg(b, "d_ips"))
    summ["d_ips2026_w_all"] = float(wavg(b, "d_ips2026"))
    summ["share_gt1km_2026"] = float(np.average(b.d_ips2026 > 1000, weights=b.persons))
    summ["d_hosp_w_all"] = float(wavg(b, "d_hosp"))
    summ["share_gt1km_all"] = float(np.average(b.d_ips > 1000, weights=b.persons))
    # 2. Moran and bivariate LISA
    lisas = {}
    for k in (8, 4, 12):
        nb = knn_weights(P, k)
        summ[f"moran_ipm_k{k}"] = moran(b.ipm.values, nb)
        summ[f"moran_dips_k{k}"] = moran(b.d_ips.values, nb)
        summ[f"bimoran_k{k}"] = global_bivariate(b.ipm.values, b.d_ips.values, nb)
        if k == 8:
            summ["p_moran_ipm_k8"] = perm_p_global(lambda x, y, n: moran(y, n), b.ipm.values, b.ipm.values, nb)
            summ["p_moran_dips_k8"] = perm_p_global(lambda x, y, n: moran(y, n), b.d_ips.values, b.d_ips.values, nb)
            summ["p_bimoran_k8"] = perm_p_global(global_bivariate, b.ipm.values, b.d_ips.values, nb)
        I, p, cl = bivariate_lisa(b.ipm.values, b.d_ips.values, nb)
        if k == 8:
            fdr = bh(p)
            lab_fdr8 = np.where(fdr, cl, "Not significant")
            summ["lisa_k8_fdr"] = {c: {"blocks": int(((cl == c) & fdr).sum()), "persons": int(b.persons[(cl == c) & fdr].sum())}
                                   for c in ["High-High", "High-Low", "Low-High", "Low-Low"]}
            summ["lisa_k8_hh_own"] = {"ipm_above_mean": int(((cl == "High-High") & (p < 0.05)).sum()),
                                      "own_distance_above_mean": int(((cl == "High-High") & (p < 0.05) & (b.d_ips > b.d_ips.mean())).sum())}
        sig = p < 0.05
        lab = np.where(sig, cl, "Not significant")
        lisas[k] = lab
        summ[f"lisa_k{k}"] = {c: {"blocks": int((lab == c).sum()), "persons": int(b.persons[lab == c].sum())}
                              for c in ["High-High", "High-Low", "Low-High", "Low-Low", "Not significant"]}
    nbh = knn_weights(P, 8)
    I, p, cl = bivariate_lisa(b.ipm.values, b.d_hosp.values, nbh)
    labh = np.where(p < 0.05, cl, "Not significant")
    summ["lisa_hosp_k8"] = {c: {"blocks": int((labh == c).sum()), "persons": int(b.persons[labh == c].sum())}
                            for c in ["High-High", "High-Low", "Low-High", "Low-Low", "Not significant"]}
    summ["bimoran_hosp_k8"] = global_bivariate(b.ipm.values, b.d_hosp.values, nbh)
    out = b[["block", "x", "y", "ipm", "persons", "d_ips", "d_hosp"]].copy()
    out["lisa_k8"] = lisas[8]
    out["lisa_k8_fdr"] = lab_fdr8
    out["lisa_hosp_k8"] = labh
    out.to_csv(RES / "lisa.csv", index=False)
    hh = out[out.lisa_k8 == "High-High"]
    summ["hh_mean_ipm"] = float(np.average(hh.ipm, weights=hh.persons)) if len(hh) else None
    summ["hh_mean_dips"] = float(np.average(hh.d_ips, weights=hh.persons)) if len(hh) else None
    # 3. formal care-seeking (main distances: sites open by mid-2018, 2018 roads)
    rows, tert = [], None
    specs = [("main", "d_ips", b, "ipm", False), ("hospital", "d_hosp", b, "ipm", False),
             ("euclidean", "e_ips", b, "ipm", False), ("snap_le_200", "d_ips", b[b.snap_m <= 200], "ipm", False),
             ("stratum", "d_ips", b, "stratum", False), ("stratum_hospital", "d_hosp", b, "stratum", False),
             ("stratum_euclidean", "e_ips", b, "stratum", False),
             ("stratum_snap_le_200", "d_ips", b[b.snap_m <= 200], "stratum", False),
             ("stratum_gp2018", "d_gp2018_osm2018", b, "stratum", False),
             ("stratum_er2018", "d_er2018_osm2018", b, "stratum", False),
             ("stratum_undated", "d_all2018_undated_osm2018", b, "stratum", False),
             ("stratum_roads2026", "d_all2018_osm2026", b, "stratum", False),
             ("stratum_network2026", "d_ips2026", b, "stratum", False),
             ("stratum_age5", "d_ips", b, "stratum", True)]
    specs += [("stratum_cluster_section", "d_ips", b, "stratum", "sec"), ("stratum_cluster_sector", "d_ips", b, "stratum", "set")]
    specs += [(f"stratum_{c[2:-8]}", c, b, "stratum", False) for c in b.columns if c.startswith("d_jit")]
    for name, dcol, sub, adj, a5 in specs:
        cut = {"sec": 2, "set": 4}.get(a5, 0) if isinstance(a5, str) else 0
        est, n, d = care_model(ill, sub, dcol, adj, age5=(a5 is True), cluster_cut=cut)
        m, X = est.attrs["model"], est.attrs["X"]
        est = pd.DataFrame(est)
        est["model"] = name
        est["n"] = n
        est["blocks"] = int(d.block.nunique())
        est["clusters"] = int(d.block.str[:len(d.block.iloc[0]) - cut].nunique()) if cut else int(d.block.nunique())
        rows.append(est.reset_index().rename(columns={"index": "term"}))
        if name == "stratum":
            # average predicted probability of formal care with the distance fixed at 0.5, 1 and 2 km
            V = m.cov_params().values
            pred = {}
            for km in (0.5, 1.0, 2.0):
                Xk = X.copy()
                Xk["dist_km"] = km
                pk = 1 / (1 + np.exp(-(Xk.values @ m.params.values)))
                grad = ((pk * (1 - pk))[:, None] * Xk.values).mean(axis=0)
                se = float(np.sqrt(grad @ V @ grad))
                pred[str(km)] = {"p": float(pk.mean()), "lo": float(pk.mean() - 1.96 * se), "hi": float(pk.mean() + 1.96 * se)}
            summ["predicted_stratum"] = pred
            # spatial autocorrelation of the residuals, aggregated by block (Pearson residuals, block means)
            r = pd.Series(m.resid_pearson, index=d.index).groupby(d.block).mean()
            br = b.set_index("block").loc[r.index]
            nbr = knn_weights(br[["x", "y"]].values, 8)
            summ["residual_moran_k8"] = moran(r.values, nbr)
            summ["p_residual_moran_k8"] = perm_p_global(lambda x, y, nn: moran(y, nn), r.values, r.values, nbr)
            summ["residual_blocks"] = int(len(r))
        if name == "main":
            # common distance tertiles for all the population
            d["tert"] = pd.qcut(d.dist_km, 3, labels=["T1 (nearest)", "T2", "T3 (farthest)"])
            tert = d.groupby("tert", observed=True).agg(persons=("formal", "size"), formal=("formal", "mean"),
                                                       dist_min=("dist_km", "min"), dist_max=("dist_km", "max"))
            iq = pd.qcut(d.ipm.rank(method="first"), 3, labels=["low", "mid", "high"])
            summ["formal_by_ipm_tertile"] = d.groupby(iq, observed=True).formal.mean().round(4).to_dict()
            dd = d.assign(iq=iq)
            rows_se = []
            for (iq_, t_), gg in dd.groupby(["iq", "tert"], observed=True):
                pbar = gg.formal.mean()
                cs = gg.groupby("block").formal.agg(["sum", "size"])
                ncl = len(cs)
                var = (ncl / (ncl - 1)) * ((cs["sum"] - pbar * cs["size"]) ** 2).sum() / len(gg) ** 2
                rows_se.append({"iq": iq_, "tert": t_, "mean": pbar, "size": len(gg), "blocks": ncl, "se": float(np.sqrt(var))})
            g = pd.DataFrame(rows_se).set_index(["iq", "tert"])  # SE clustered by block
            g.reset_index().to_csv(RES / "care_ipm_by_distance.csv", index=False)
            summ["ill_in_model"] = int(n)
            summ["formal_share_model"] = float(d.formal.mean())
    jit = pd.concat(rows)
    jit = jit[jit.model.str.startswith("stratum_jit") & (jit.term == "dist_km")]
    summ["jitter_or_range"] = {r: [float(jit[jit.model.str.contains(f"jit{r}_")]["or"].min()),
                                   float(jit[jit.model.str.contains(f"jit{r}_")]["or"].max())] for r in (100, 250)}
    summ["jitter_or_p2575"] = {r: [float(jit[jit.model.str.contains(f"jit{r}_")]["or"].quantile(q_)) for q_ in (0.025, 0.5, 0.975)]
                               for r in (100, 250)}
    summ["jitter_draws"] = int(jit.model.str.contains("jit100_").sum())
    summ["jitter_spearman_range"] = {r: [float(min(spearmanr(b.ipm, b[c])[0] for c in b.columns if c.startswith(f"d_jit{r}_"))),
                                         float(max(spearmanr(b.ipm, b[c])[0] for c in b.columns if c.startswith(f"d_jit{r}_")))]
                                     for r in (100, 250)}
    pd.concat(rows).to_csv(RES / "care_models.csv", index=False)
    tert.to_csv(RES / "care_tertiles.csv")
    (RES / "summary.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
    print(q.round(3).to_string())
    print(json.dumps({k: v for k, v in summ.items() if not k.startswith("lisa_k4") and not k.startswith("lisa_k12")},
                     ensure_ascii=False, indent=1))
    print(pd.concat(rows).round(4).to_string())
    print(tert.round(4))


if __name__ == "__main__":
    main()
