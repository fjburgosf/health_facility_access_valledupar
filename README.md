# Pobreza multidimensional, logro educativo, distancia a los servicios de salud y búsqueda de atención en las manzanas de Valledupar, Colombia

Este repositorio contiene el código de análisis del artículo sometido a *Territorios*. En el estudio se examinó si la pobreza multidimensional, el bajo logro educativo adulto y la distancia a los servicios de salud se concentran en las mismas manzanas urbanas de Valledupar (Cesar, Colombia). También se analizó si la distancia se asocia con la búsqueda de atención formal entre las personas que declararon un problema de salud en el Censo Nacional de Población y Vivienda de 2018 (CNPV).

Para estimar la distancia se utilizaron las sedes de instituciones prestadoras de servicios de salud (IPS) con fecha de apertura hasta el 1 de julio de 2018, según las fechas que figuran en el Registro Especial de Prestadores de Servicios de Salud (REPS) consultado en 2026. Se utilizó la red vial de OpenStreetMap del 1 de julio de 2018. Las distancias con la red y las sedes de 2026 se calcularon para la comparación.

Se incluyen los scripts de análisis, dos archivos con datos institucionales del REPS, el registro de procedencia de las fuentes descargadas ([`data/raw/PROVENANCE.json`](data/raw/PROVENANCE.json), con sus valores SHA-256) y los [resultados agregados](results/). Los microdatos censales y los resultados por manzana se utilizan durante la ejecución local.

## Datos

| Fuente | Obtención | Ubicación |
|---|---|---|
| Índice de pobreza multidimensional (IPM) censal por manzana, DANE (2018) | Descarga en el paso 01 | `data/raw/ipm_manzana_valledupar.geojson` |
| Equipamientos de salud y vías de OpenStreetMap (consulta actual) | Descarga en el paso 01 mediante la API de Overpass | `data/raw/osm_health.json`, `data/raw/osm_roads.json` |
| Red vial de OpenStreetMap al 1 de julio de 2018 | Consulta histórica en el paso 01b mediante la API de Overpass | `data/raw/osm_roads_2018.json` |
| Sedes del REPS (datos abiertos, conjunto `c36g-9fc2` de datos.gov.co) | Descarga en el paso 01. Se excluyen los registros de profesionales independientes | `data/raw/reps_valledupar_institutional.json` |
| Fechas de apertura y servicios habilitados de las sedes de IPS de Valledupar | [Portal del REPS](https://prestadores.minsalud.gov.co/habilitacion/), módulo «Registro actual», con acceso de invitado. Se incluyen los campos institucionales consultados el 21 de noviembre de 2026 | `data/raw/reps_portal_2026-11-21/` |
| Microdatos del CNPV 2018 para el Cesar | Catálogo del Archivo Nacional de Datos (ANDA) del DANE, previo registro | `data/cnpv2018_cesar_csv/` (`CNPV2018_5PER_A2_20.CSV`, `CNPV2018_1VIV_A2_20.CSV`, `CNPV2018_MGN_A2_20.CSV`) |

En `servicios_flags_sede.txt` se codifican los servicios habilitados de cada sede mediante una máscara de bits. El valor 1 corresponde a consulta externa de medicina general o familiar (servicios 328 o 325), el 2 a urgencias (1102) y el 4 a hospitalización de adultos o pediátrica (129 o 130).

## Instalación

El análisis se ejecutó con Python 3.12.

```bash
python -m venv .venv
# Activación en Windows: .venv\Scripts\activate
# Activación en Linux o macOS: source .venv/bin/activate
pip install -r requirements.txt
```

## Ejecución del análisis

Ejecute los scripts desde la raíz del repositorio en el orden indicado.

```bash
python src/01_get_sources.py
python src/01b_get_roads_2018.py
python src/02_geocode_facilities.py
python src/03_blocks_network.py
python src/03b_network_2018.py
python src/04_analysis.py
python src/04b_education.py
```

El paso 01 descarga las versiones disponibles al momento de la consulta y el paso 01b recupera la red vial histórica. Los archivos obtenidos en una descarga posterior podrían diferir de los utilizados en el estudio. Cada ejecución registra los valores SHA-256 de los archivos descargados en `data/raw/PROVENANCE.json`.

| Paso | Script | Operación | Salida principal |
|---|---|---|---|
| 01 | `01_get_sources.py` | Descarga del IPM por manzana, los equipamientos y las vías de OpenStreetMap y los datos abiertos del REPS | `data/raw/` |
| 01b | `01b_get_roads_2018.py` | Descarga de la red vial de OpenStreetMap al 1 de julio de 2018 | `data/raw/osm_roads_2018.json` |
| 02 | `02_geocode_facilities.py` | Geocodificación local de las direcciones de las sedes de IPS a partir de la nomenclatura vial | `data/processed/facilities.csv` |
| 03 | `03_blocks_network.py` | Integración de las manzanas urbanas con la población y las respuestas censales de salud. Construcción del grafo vial y cálculo de distancias para 2026 | `data/processed/blocks.csv`, `data/processed/ill_persons.csv` |
| 03b | `03b_network_2018.py` | Reconstrucción de la red de 2018 (todas las sedes, consulta de medicina general, urgencias y hospitales o clínicas) y desplazamiento aleatorio de las sedes | `data/processed/blocks_2018.csv` |
| 04 | `04_analysis.py` | Comparación de distancias por quintil del IPM, correlaciones de Spearman, estadísticos de Moran global y local bivariado, modelos logísticos de búsqueda de atención y análisis de sensibilidad | `results/summary.json`, `results/care_models.csv`, `results/quintiles.csv`, `results/lisa.csv` |
| 04b | `04b_education.py` | Análisis del bajo logro educativo adulto por manzana y del clima educativo del hogar, su concentración espacial conjunta con la distancia y los modelos de búsqueda de atención con interacción | `results/education_*.csv`, `results/education_summary.json` |

## Definiciones

- **Pobreza por manzana.** IPM censal de la manzana (porcentaje de personas que viven en hogares en situación de pobreza multidimensional).
- **Atención formal.** Entre quienes declararon un problema de salud en los 30 días previos al censo (`P_ENFERMO` = 1), se consideró atención formal acudir a la entidad aseguradora de salud o a un médico particular (`P_QUEHIZO_PPAL` = 1 o 2). Los códigos 3 a 9 se agruparon como otras respuestas y se excluyeron los códigos 0, 11 y 99.
- **Estrato.** Estrato socioeconómico de la vivienda (`VA1_ESTRATO`).
- **Logro educativo.** Los niveles de `P_NIVEL_ANOSR` se agruparon en primaria o menos (10, 1 y 2), secundaria (3 a 6) y educación superior (7 a 9). El clima educativo del hogar se definió como el nivel más alto alcanzado por sus miembros de 25 años o más.
- **Distancia por la red vial.** Se calculó con el algoritmo de Dijkstra de múltiples orígenes sobre el mayor componente conexo del grafo vial no dirigido. Se sumaron las distancias de conexión de las manzanas y de las sedes a los nodos de la red.

## Análisis de sensibilidad e incertidumbre

Los estadísticos locales de Moran se evaluaron con 999 permutaciones condicionales (semilla 20261120). Los vecinos de cada permutación se seleccionaron sin reemplazo y la tasa de falsos descubrimientos se controló con el procedimiento de Benjamini-Hochberg. Si varias sedes se conectaban al mismo nodo vial, se conservó la mínima distancia de conexión. En el análisis de sensibilidad con desplazamientos aleatorios de las sedes se generaron 50 realizaciones por radio. Los errores estándar de los modelos logísticos se agruparon por manzana y, en los análisis de sensibilidad, por sección y por sector censal. Las probabilidades predichas se estimaron con intervalos de confianza obtenidos mediante el método delta.
