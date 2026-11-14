"""AnnotateMissense: lookup, VCF annotation and serving of precomputed
missense pathogenicity predictions for GRCh38/hg38.

The core package depends only on DuckDB. The optional web server lives in the
top-level ``app`` package and additionally requires Flask.
"""

__version__ = "1.1.0"

# Canonical names used by the API, CLI and VCF output. The released database
# keeps its historical column names (UQ_score, UQ_sd, UQ_pred); see
# docs/naming.md in the AnnotateMissense repository.
SCORE_FIELD = "annotatemissense_score"
CONFIDENCE_FIELD = "annotatemissense_confidence"
PRED_FIELD = "annotatemissense_pred"

LEGACY_COLUMN_NAMES = {
    "UQ_score": SCORE_FIELD,
    "UQ_sd": CONFIDENCE_FIELD,
    "UQ_pred": PRED_FIELD,
}

DECISION_THRESHOLD = 0.5
GENOME_BUILD = "GRCh38"
ZENODO_DOI = "10.5281/zenodo.19981867"
