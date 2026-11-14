# data/

Put the AnnotateMissense data files here (they are not part of the git repository).

| File | Source | Used for |
|---|---|---|
| `variants.duckdb` | `variants.duckdb.gz` from <https://doi.org/10.5281/zenodo.19981867> (`gzip -dk variants.duckdb.gz`) | full annotations + scores (337 columns) |
| `annotatemissense_scores.duckdb` | built from `UQ_BioSig_model_Final.tsv.gz` with `annotatemissense build-score-index` | genome-wide scores only (light option) |
| `variants_dev.duckdb` | chromosome 22 extract of `variants.duckdb` | demo / development |

The server picks `variants.duckdb`, then `variants_dev.duckdb`, then `annotatemissense_scores.duckdb`
from this folder automatically, or use the `ANNOTATEMISSENSE_DB` / `ANNOTATEMISSENSE_SCORES`
environment variables.
