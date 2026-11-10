# Custody architecture of the public wastewater sequencing archive

Record-level data, analysis code and results for the audit of submission
architecture (ownership, routing and submission velocity) across the three
strata of the public INSDC wastewater record: the SARS-CoV-2 census, the
non-SARS pathogen census and the wastewater-metagenome random sample.

## Contents

- `data/` — audited record-level tables (CC0):
  - `master_biosample.csv` — SARS-CoV-2 wastewater census (82,757 records,
    NCBI BioSample, ascertained 15 September 2026)
  - `multipathogen_biosample.csv` — non-SARS pathogen census (1,426 records;
    influenza A, influenza B, norovirus, rotavirus A, RSV; ascertained
    2 November 2026)
  - `metagenome_sample.csv` — metagenome random sample (3,000 of 373,151
    records; ascertained 2 November 2026)
  - `who_covid_global.csv` — WHO COVID-19 dashboard extract (23 August 2026)
  - `country_stats.csv`, `burden_country.csv`, `misalignment_matrix.csv` —
    derived per-country tables
- `code/` — acquisition (`fetch_full.py`, `fetch_multipathogen.py`),
  country attribution and misalignment analysis (`map_analyze.py`,
  `analyze_multipathogen.py`), statistical models (`upgrade_stats.py`),
  reviewer-response reanalyses (`reanalysis.py`, `bootstrap_fast.py`,
  `aft_fit.py`) and figure code (`make_flowfig.py`)
- `results/` — JSON outputs of all models (Cox refits, AFT, RMST, Firth and
  drop-Peru presence models, cluster bootstrap, record accounting)
- `REQUIREMENTS.txt` — pinned Python and library versions
- `Supplementary_Data_Dictionary.md` — data dictionary for all tables

## Reproducibility

- The metagenome sample was drawn with Python `random.sample`, fixed seed 42
  (`code/fetch_multipathogen.py`).
- Census query: `SARS-CoV-2[Organism] AND wastewater[All Fields]` (NCBI
  BioSample, E-utilities). Synonym sensitivity queries (sewage, sewer, sludge;
  influent; effluent) are documented in the manuscript Methods.
- All analyses run on Python 3.13 with the versions pinned in
  `REQUIREMENTS.txt`.
- Record accounting (Supplementary Figure 3 of the manuscript) is reproduced
  by `code/reanalysis.py` (section "Part 1").

## License

Record-level tables and derived data: CC0. Code: MIT.
