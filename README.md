# Poverty, education and distance to health facilities in the city blocks of Valledupar, Colombia (2018): analysis code

This repository contains the analysis code used to study the joint spatial distribution of multidimensional poverty,
adult educational attainment and network distance to health facilities in the urban blocks of Valledupar (Cesar,
Colombia), and the association between distance and formal care-seeking among persons with a health problem in the
2018 National Population and Housing Census (CNPV).

Distances are computed with the health facilities that were open at the time of the census (sites of health care
institutions, IPS, open by 1 July 2018 according to the national registry of health providers, REPS) and with the
OpenStreetMap road network as it was on that date. The 2026 network is kept for comparison.

The repository includes only the analysis scripts and two small files with institutional data of the REPS. It does
not include the census microdata.

## Data

| Source | How to obtain it | Location |
|---|---|---|
| IPM by city block, DANE (census multidimensional poverty, 2018) | downloaded by step 01 | `data/raw/ipm_manzana_valledupar.geojson` |
| OpenStreetMap health facilities and roads (current) | downloaded by step 01 (Overpass API) | `data/raw/osm_health.json`, `osm_roads.json` |
| OpenStreetMap roads on 1 July 2018 | downloaded by step 01b (Overpass API, historical query) | `data/raw/osm_roads_2018.json` |
| REPS sites (open data, datos.gov.co `c36g-9fc2`) | downloaded by step 01; independent professionals are not stored | `data/raw/reps_valledupar_institutional.json` |
| REPS query module: opening date and enabled services of the IPS sites of Valledupar | REPS portal (<https://prestadores.minsalud.gov.co/habilitacion/>, "Registro actual", guest access), retrieved on 2026-11-21; included here (institutional fields only) | `data/raw/reps_portal_2026-11-21/` |
| CNPV 2018 microdata, Cesar | DANE microdata catalogue (ANDA), after registration; not redistributed | `data/cnpv2018_cesar_csv/` (`CNPV2018_5PER_A2_20.CSV`, `CNPV2018_1VIV_A2_20.CSV`, `CNPV2018_MGN_A2_20.CSV`) |

`servicios_flags_sede.txt` summarises the enabled services of each site as a bit mask: 1 = general or family
medicine outpatient consultation (services 328 or 325), 2 = emergency service (1102), 4 = adult or paediatric
hospitalisation (129, 130).

## Installation

Python 3.12 was used.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate    Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

## Running the analysis

Run the scripts from the repository root, in this order:

```bash
python src/01_get_sources.py
python src/01b_get_roads_2018.py
python src/02_geocode_facilities.py
python src/03_blocks_network.py
python src/03b_network_2018.py
python src/04_analysis.py
python src/04b_education.py
```

Steps 01 and 01b download current versions of online sources, so a later download could differ slightly from the
files used in the study (their SHA-256 values are recorded in `data/raw/PROVENANCE.json`).

| Step | Script | What it does | Main output |
|---|---|---|---|
| 01 | `01_get_sources.py` | Downloads the block IPM, OpenStreetMap facilities and roads, and the REPS open data | `data/raw/` |
| 01b | `01b_get_roads_2018.py` | Downloads the OpenStreetMap roads as they were on 1 July 2018 | `data/raw/osm_roads_2018.json` |
| 02 | `02_geocode_facilities.py` | Local geocoding of the IPS addresses on the street grid (no external geocoder) | `data/processed/facilities.csv` |
| 03 | `03_blocks_network.py` | Urban blocks, census population and health questions, road graph and network distances (2026) | `data/processed/blocks.csv`, `ill_persons.csv` |
| 03b | `03b_network_2018.py` | Network of 2018 (sites open by mid-2018, general medicine, emergency, hospitals), 2018 roads, random displacement of the sites | `data/processed/blocks_2018.csv` |
| 04 | `04_analysis.py` | Distance by IPM quintile, Spearman correlations, global and bivariate local Moran (with false discovery rate), logistic models of formal care-seeking and sensitivity analyses | `results/summary.json`, `care_models.csv`, `quintiles.csv`, `lisa.csv` |
| 04b | `04b_education.py` | Low adult education by block, household educational climate, co-accumulation and care-seeking models with interaction | `results/education_*.csv`, `education_summary.json` |

## Definitions

- Block poverty: census IPM of the block (percentage of persons in multidimensionally poor households).
- Formal care: persons with a health problem in the last 30 days (`P_ENFERMO` = 1) who went to their health insurer
  or to a private doctor (`P_QUEHIZO_PPAL` in 1–2); responses 3–9 are other responses; 0, 11 and 99 are excluded.
- Stratum: socioeconomic stratum of the dwelling (`VA1_ESTRATO`).
- Education (`P_NIVEL_ANOSR`): primary or less (10, 1, 2), secondary (3–6), tertiary (7–9). Household educational
  climate: highest level among the members aged 25 or more.
- Network distance: multi-source Dijkstra on the undirected road graph (largest connected component), plus the
  snapping distances of blocks and sites.

## Uncertainty

Local Moran statistics use 999 conditional permutations (seed 20261120) and the Benjamini-Hochberg false discovery
rate. Logistic models use standard errors clustered by block (and by census section and sector in sensitivity
analyses). Predicted probabilities use the delta method.
