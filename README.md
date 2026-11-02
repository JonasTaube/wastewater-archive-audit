# wastewater-archive-audit

Code accompanying the Article *"Same sewers, different sovereignty: a stratified audit of
457,334 public wastewater sequencing records shows submission architecture - not capacity,
income, or disease burden - shapes the global archive"* (Tianxi Yu et al., 2026; v2 of the
single-gateway audit).

## What the audit measures

Three strata of the public INSDC wastewater record:

| Stratum | Records | Method | EBI gateway | Domestic accounts | Countries |
|---|---|---|---|---|---|
| SARS-CoV-2 wastewater | 82,757 | complete census (2026-09-15) | 98.9% | 1.1% | 20 |
| Other pathogen targets (influenza A/B, norovirus, rotavirus A, RSV) | 1,426 | complete census (2026-11-02) | 1.2% | 76.5% | 9 |
| Wastewater metagenome | 373,151 | uniform random sample, n=3,000 (2026-11-02) | 19.8% | 60.7% | 57 |

Key statistical upgrades in v2: Cox proportional-hazards model of submission velocity
(hub-routed HR 0.093, 95% CI 0.085-0.101) and a two-part presence model over 213 WHO
entities (UMIC presence OR 0.032 vs HIC at equal mortality and population).

## Pipeline

| Script | Purpose |
|---|---|
| `fetch_pilot.py` | Pilot retrieval: query formulation and union sensitivity query (wastewater/sewage/sewer/sludge). |
| `fetch_full.py` | Full retrieval of all SARS-CoV-2 records (accession, Owner, submission date, geo location, collection date, organism, isolation source, BioProject). |
| `fetch_multipathogen.py` | **v2:** full retrieval of the non-SARS pathogen census (5 organisms) and the 3,000-record random sample of the 373,151-record wastewater metagenome stratum. |
| `map_analyze.py` | Country attribution (three-tier rule system) and core descriptive statistics for the SARS-CoV-2 stratum. |
| `analyze_multipathogen.py` | **v2:** stratified analysis of the non-SARS pathogen and metagenome strata (extended institution dictionary: NWSS, state public-health laboratories, Guangdong CDC, Aalborg, etc.). |
| `burden_analysis.py` | Linkage to the WHO COVID-19 Dashboard cumulative mortality series. |
| `reconcile_and_burden.py` | Flow reconciliation across pools; two-level burden analysis; back-filling sensitivity. |
| `upgrade_stats.py` | **v2:** Cox proportional-hazards velocity model (lifelines) and two-part presence model (statsmodels, World Bank covariates). |
| `owner_semantics_check.py` | Verification of BioSample Owner semantics against the ENA submission registry. |
| `make_figures_nm.py` | Figure rendering (matplotlib, Arial, 300 dpi, PDF fonttype 42). |
| `resolve_nm_refs.py` | Reference list resolution and citation-continuity checks. |

## Results (JSON)

- `multipathogen_fetch.json` - stratum ascertainment counts and sampling design.
- `multipathogen_analysis.json` - stratified custody/lag/income analysis for strata 2 and 3.
- `stats_upgrade.json` - Cox model and two-part presence model estimates.

## Requirements

Python 3.13; `pandas`, `numpy`, `scipy`, `matplotlib`, `requests`, `statsmodels`, `lifelines`.
`pip install -r requirements.txt`

## Data (included in this repository)

All audited record-level tables are openly available under [`data/`](data/):

| File | Content | Records |
|---|---|---|
| `data/master_biosample.csv` | SARS-CoV-2 wastewater census with country attributions | 82,757 |
| `data/multipathogen_biosample.csv` | Non-SARS pathogen census (influenza A/B, norovirus, rotavirus A, RSV) | 1,426 |
| `data/metagenome_sample.csv` | Wastewater-metagenome random sample | 3,000 of 373,151 |
| `data/who_covid_global.csv` | WHO COVID-19 mortality series (burden linkage) | 239 entities |
| `data/burden_country.csv` | Country-level burden linkage | 20 contributors |
| `data/country_stats.csv` / `data/misalignment_matrix.csv` | Attribution statistics / flow matrix | — |

Verification results are under [`results/`](results/). Source records are publicly accessible
through NCBI BioSample; a Zenodo DOI archive is minted on acceptance.

## Licence

MIT - see [LICENSE](LICENSE).
