# wastewater-archive-audit

Code accompanying the Article *"A single national gateway: a global audit of the public
wastewater sequencing archive reveals an inequitable submission architecture misaligned
with disease burden"* (Tianxi Yu, 2026).

The audit is a complete census of all 82,757 SARS-CoV-2 wastewater sequencing records in
the NCBI BioSample archive (ascertained 15 September 2026), attributing every record to
its submitting institution and submission pathway, and linking the result to WHO COVID-19
mortality data.

## Pipeline

| Script | Purpose |
|---|---|
| `fetch_pilot.py` | Pilot retrieval: NCBI E-utilities esearch/efetch query formulation and union sensitivity query (wastewater/sewage/sewer/sludge). |
| `fetch_full.py` | Full retrieval of all records (accession, Owner, submission date, geo location, collection date, organism, isolation source, BioProject). |
| `map_analyze.py` | Country attribution (three-tier rule system: transnational database proxies; curated institution dictionary; city-to-country dictionary) and core descriptive statistics. |
| `burden_analysis.py` | Linkage to the WHO COVID-19 Dashboard cumulative mortality series; burden-normalised intensity and representation-mismatch ratios. |
| `reconcile_and_burden.py` | Flow reconciliation of all 82,757 records across pools; two-level burden analysis (all 239 WHO entities and 20 contributing countries); back-filling sensitivity analyses. |
| `owner_semantics_check.py` | Verification of BioSample Owner semantics: 300-record random sample checked against the ENA submission registry (centre names). |
| `make_figures_nm.py` | Figure rendering (matplotlib, Arial, 300 dpi, PDF fonttype 42). |
| `resolve_nm_refs.py` | Reference list resolution and citation-continuity checks. |

## Requirements

Python 3.13; `pandas`, `numpy`, `scipy`, `matplotlib`, `requests`.
`pip install -r requirements.txt`

## Data

The full audited record-level table (Supplementary Data 1) is deposited on Zenodo
(DOI minted on acceptance). Source records are publicly accessible through NCBI
BioSample under the BioProject accessions listed in that table.

## Licence

MIT — see [LICENSE](LICENSE).
