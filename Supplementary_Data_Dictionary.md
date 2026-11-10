# Data Dictionary — Supplementary Data 1–3

Wastewater archive audit (Sci Data submission). All fields derive from public
INSDC (NCBI BioSample / ENA) metadata as of the ascertainment dates given in
Methods, plus the WHO COVID-19 dashboard extract (23 August 2026).

## Supplementary Data 1 — SARS-CoV-2 census (82,757 records)

| Field | Type | Description |
|---|---|---|
| `accession` | string | BioSample accession (SAMN* = NCBI-native; SAMEA* = ENA-native, reached BioSample via the INSDC partnership). |
| `owner` | string | BioSample Owner element: the submitting account / institution of record. |
| `submission_date` | ISO 8601 (may be truncated to month) | Date the record was made public in the archive. |
| `geo_location` | string | `geo_loc_name` attribute, `Country: region` format. |
| `collection_date` | ISO 8601 (may be truncated to month) | Sample collection date. |
| `organism` | string | NCBI taxonomy name (constant within each stratum). |
| `isolation_source` | string | Submitter-provided source description. |
| `project` | string | BioProject accession. |
| `sample_title` | string | Submitter-provided title. |

Derived fields (added by the audit pipeline):
| `sample_country` | ISO 3166-1 alpha-3 | First token of `geo_location`, harmonised; `(missing)` for the 84 records without geolocation. |
| `owner_country` | ISO 3166-1 alpha-3 | Country of registration of the submitting institution (three-tier rule system: database proxies → curated institution dictionary → city dictionary). |

## Supplementary Data 2 — non-SARS pathogen census (1,426 records)

Same schema as Data 1; `organism` ∈ {Influenza A virus, Influenza B virus,
Norovirus, Rotavirus A, Respiratory syncytial virus}; `pathogen_layer` flags
the stratum.

## Supplementary Data 3 — metagenome random sample (3,000 records)

Same schema as Data 1; drawn by uniform random sampling without replacement
from the 373,151 wastewater-metagenome records (`wastewater metagenome[Organism]`,
ascertained 2 November 2026) with fixed seed 42 (Python `random.sample`,
`fetch_multipathogen.py`); sampling fraction 0.8%.

## Analysis outputs

| File | Content |
|---|---|
| `country_stats.csv` | Per-sampling-country statistics (counts, misalignment, foreign-submission rate, median lag). |
| `misalignment_matrix.csv` | Sampling country × submitting-institution country cross-tabulation (foreign-attributed records). |
| `robustness_stats.json` | Cox model, Schoenfeld tests, clustered SEs, leave-one-country-out results. |
| `reviewer_response_stats.json` | Revised analyses: window exclusions, no-window and complete-follow-up Cox refits, AFT, RMST, Firth/drop-Peru presence models, cluster bootstrap. |
| `burden_country.csv` | Country-level archive share vs WHO cumulative deaths (Table 1 source). |
