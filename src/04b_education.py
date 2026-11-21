"""Step 04b: educational attainment, distance and formal care-seeking.

Distances: IPS sites open by mid-2018 on the 2018 road network, as in step 04.

Definitions (CNPV 2018, P_NIVEL_ANOSR; codebook cnpv2018.ddi.xml):
- low attainment = none (10), preschool (1) or basic primary (2); secondary = basic secondary, upper secondary
  (academic or technical) or normal school (3-6); tertiary = technical, technological, university or postgraduate
  (7-9). "No aplica" (12), "no informa" (99) and missing are excluded.
- Block low adult education: share of the residents aged 25 or more with low attainment.
- Household educational climate: highest attainment among the household members aged 25 or more
  (household = COD_ENCUESTAS + U_VIVIENDA + P_NROHOG). It applies to every member, including children.

Analyses:
1. Spearman correlation between block low adult education and network distance; population-weighted distance by
   quintile of low adult education; global and local bivariate Moran (k = 8, 999 permutations, p < 0.05), with
   the same functions as step 04.
2. Formal care-seeking: the stratum-adjusted model of step 04 plus household educational climate (reference:
   tertiary), and the same model with a distance x climate interaction, from which the OR per km within each
   climate group is obtained.

Outputs (results): education_summary.json, education_quintiles.csv, education_care_models.csv,
education_lisa.csv

Usage (from the repository root):  python src/04b_education.py
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DAT = ROOT / "data" / "processed"
RES = ROOT / "results"
CNPV = ROOT / "data" / "cnpv2018_cesar_csv"
KEY = ["COD_ENCUESTAS", "U_VIVIENDA", "P_NROHOG"]
LEVEL = {"10": 0, "1": 0, "2": 0, "3": 1, "4": 1, "5": 1, "6": 1, "7": 2, "8": 2, "9": 2}
CLIMATE = {0: "primary or less", 1: "secondary", 2: "tertiary"}

_spec = importlib.util.spec_from_file_location("analysis04", Path(__file__).with_name("04_analysis.py"))
a04 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a04)


def census():
    mgn = pd.read_csv(CNPV / "CNPV2018_MGN_A2_20.CSV", dtype=str,
                      usecols=["U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA", "COD_DANE_ANM"])
    mgn = mgn[(mgn.U_MPIO == "001") & (mgn.UA_CLASE == "1")]
    per = pd.read_csv(CNPV / "CNPV2018_5PER_A2_20.CSV", dtype=str,
                      usecols=["U_MPIO", "UA_CLASE", *KEY, "P_SEXO", "P_EDADR", "P_ENFERMO", "P_QUEHIZO_PPAL",
                               "P_NIVEL_ANOSR"])
    per = per[(per.U_MPIO == "001") & (per.UA_CLASE == "1")]
    per = per.merge(mgn[["COD_ENCUESTAS", "U_VIVIENDA", "COD_DANE_ANM"]], on=["COD_ENCUESTAS", "U_VIVIENDA"], how="left")
    viv = pd.read_csv(CNPV / "CNPV2018_1VIV_A2_20.CSV", dtype=str,
                      usecols=["U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA", "VA1_ESTRATO"])
    viv = viv[(viv.U_MPIO == "001") & (viv.UA_CLASE == "1")]
    per = per.merge(viv[["COD_ENCUESTAS", "U_VIVIENDA", "VA1_ESTRATO"]], on=["COD_ENCUESTAS", "U_VIVIENDA"], how="left")
    per["level"] = per.P_NIVEL_ANOSR.map(LEVEL)
    per["adult25"] = pd.to_numeric(per.P_EDADR) >= 6
    return per


def fit(d, X):
    X = sm.add_constant(X)
    groups = pd.factorize(d.block)[0]
    return sm.GLM(d.formal, X, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": groups})


def or_table(m):
    ci = m.conf_int()
    return pd.DataFrame({"or": np.exp(m.params), "lo": np.exp(ci[0]), "hi": np.exp(ci[1]), "p": m.pvalues})


def main():
    per = census()
    adults = per[per.adult25 & per.level.notna()]
    clim = adults.groupby(KEY)["level"].max().rename("climate_lvl")
    per = per.merge(clim, left_on=KEY, right_index=True, how="left")
    low = adults.assign(low=(adults.level == 0).astype(float)).groupby("COD_DANE_ANM").agg(
        low_edu=("low", "mean"), adults_valid=("low", "size"))

    b = a04.load_blocks()  # main distances: IPS sites open by mid-2018, 2018 roads
    b = b.merge(low, left_on="block", right_index=True, how="left")
    summ = {"adults25_total": int(per.adult25.sum()), "adults25_valid_level": int(len(adults)),
            "adults25_low_share": float((adults.level == 0).mean()),
            "blocks": int(len(b)), "blocks_without_adult_level": int(b.low_edu.isna().sum())}
    be = b[b.low_edu.notna()].reset_index(drop=True)
    be["low_edu_pct"] = 100 * be.low_edu
    rho, p = spearmanr(be.low_edu, be.d_ips)
    summ["spearman_lowedu_d_ips"] = {"rho": float(rho), "p": float(p)}
    rho, p = spearmanr(be.low_edu, be.ipm)
    summ["spearman_lowedu_ipm"] = {"rho": float(rho), "p": float(p)}
    wavg = lambda g, c: np.average(g[c], weights=g.persons)
    be["edu_q"] = pd.qcut(be.low_edu.rank(method="first"), 5, labels=[1, 2, 3, 4, 5])
    q = be.groupby("edu_q", observed=True).apply(lambda g: pd.Series({
        "blocks": len(g), "persons": g.persons.sum(), "low_edu_min": 100 * g.low_edu.min(),
        "low_edu_max": 100 * g.low_edu.max(), "low_edu_mean_w": 100 * wavg(g, "low_edu"),
        "d_ips_w": wavg(g, "d_ips"), "share_gt1km": np.average(g.d_ips > 1000, weights=g.persons)}),
        include_groups=False)
    q.to_csv(RES / "education_quintiles.csv")
    P = be[["x", "y"]].values
    nb = a04.knn_weights(P, 8)
    summ["moran_lowedu_k8"] = a04.moran(be.low_edu.values, nb)
    summ["bimoran_lowedu_k8"] = a04.global_bivariate(be.low_edu.values, be.d_ips.values, nb)
    I, pv, cl = a04.bivariate_lisa(be.low_edu.values, be.d_ips.values, nb)
    lab = np.where(pv < 0.05, cl, "Not significant")
    summ["lisa_lowedu_k8"] = {c: {"blocks": int((lab == c).sum()), "persons": int(be.persons[lab == c].sum())}
                              for c in ["High-High", "High-Low", "Low-High", "Low-Low", "Not significant"]}
    hh = lab == "High-High"
    summ["hh_lowedu_mean"] = float(100 * np.average(be.low_edu[hh], weights=be.persons[hh]))
    summ["hh_lowedu_dips"] = float(np.average(be.d_ips[hh], weights=be.persons[hh]))
    lisa_ipm = pd.read_csv(RES / "lisa.csv", dtype={"block": str})[["block", "lisa_k8"]]
    both = be[["block", "persons"]].assign(lisa_edu=lab).merge(lisa_ipm, on="block", how="left")
    summ["hh_overlap_with_ipm_hh"] = {"blocks": int(((both.lisa_edu == "High-High") & (both.lisa_k8 == "High-High")).sum()),
                                      "persons": int(both.persons[(both.lisa_edu == "High-High") & (both.lisa_k8 == "High-High")].sum())}
    be.assign(lisa_lowedu_k8=lab)[["block", "x", "y", "low_edu_pct", "persons", "d_ips", "lisa_lowedu_k8"]].to_csv(
        RES / "education_lisa.csv", index=False)

    # 2. care-seeking: same ill-person selection as step 03
    ill = per[(per.P_ENFERMO == "1") & per.P_QUEHIZO_PPAL.isin([str(k) for k in range(1, 10)])].copy()
    ill["formal"] = ill.P_QUEHIZO_PPAL.isin(["1", "2"]).astype(int)
    ill = ill.rename(columns={"COD_DANE_ANM": "block"})
    summ["ill_with_response"] = int(len(ill))
    _, n_ref, d = a04.care_model(ill, b, "d_ips", "stratum")
    summ["ill_in_stratum_model"] = int(n_ref)
    d = d[d.climate_lvl.notna()].copy()
    summ["ill_with_climate"] = int(len(d))
    d["climate"] = pd.Categorical(d.climate_lvl.map(CLIMATE), categories=["tertiary", "secondary", "primary or less"])
    summ["ill_by_climate"] = d.climate.value_counts().to_dict()
    summ["formal_by_climate"] = d.groupby("climate", observed=True).formal.mean().round(4).to_dict()
    summ["age65_share_by_climate"] = d.groupby("climate", observed=True).age_grp.apply(lambda s: float((s == "65+").mean())).round(4).to_dict()
    base = pd.get_dummies(d[["dist_km", "female", "age_grp", "stratum", "climate"]],
                          columns=["age_grp", "stratum", "climate"], drop_first=True, dtype=float)
    m1 = fit(d, base)
    t1 = or_table(m1).assign(model="stratum_climate", n=len(d))
    Xi = base.copy()
    for g in ("secondary", "primary or less"):
        Xi[f"dist_x_{g}"] = Xi["dist_km"] * Xi[f"climate_{g}"]
    m2 = fit(d, Xi)
    t2 = or_table(m2).assign(model="stratum_climate_interaction", n=len(d))
    # OR per km within each climate group, from linear combinations of the interaction model
    names = list(m2.params.index)
    groups = {}
    for g in ("tertiary", "secondary", "primary or less"):
        c = np.zeros(len(names))
        c[names.index("dist_km")] = 1
        if g != "tertiary":
            c[names.index(f"dist_x_{g}")] = 1
        est = float(c @ m2.params.values)
        se = float(np.sqrt(c @ m2.cov_params().values @ c))
        groups[g] = {"or": float(np.exp(est)), "lo": float(np.exp(est - 1.96 * se)), "hi": float(np.exp(est + 1.96 * se)),
                     "n": int((d.climate == g).sum())}
    summ["or_km_by_climate"] = groups
    R = np.eye(len(names))[[names.index("dist_x_secondary"), names.index("dist_x_primary or less")]]
    w = m2.wald_test(R, scalar=True)
    summ["interaction_wald"] = {"chi2": float(w.statistic), "df": 2, "p": float(w.pvalue)}
    tabs = pd.concat([t1, t2]).reset_index().rename(columns={"index": "term"})
    tabs.to_csv(RES / "education_care_models.csv", index=False)
    (RES / "education_summary.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
    print(q.round(3).to_string())
    print(json.dumps(summ, ensure_ascii=False, indent=1))
    print(tabs.round(4).to_string())


if __name__ == "__main__":
    main()
