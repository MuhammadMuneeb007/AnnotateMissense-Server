"""Read-only access to the AnnotateMissense DuckDB resources.

Two resources are supported, separately or together:

* the **annotation database** (``variants.duckdb`` from Zenodo, table
  ``variants``): 337 annotation columns plus the AnnotateMissense prediction,
  keyed by ``Chr``/``Start``/``Ref``/``Alt``;
* the **score index** (``annotatemissense_scores.duckdb``, table ``scores``):
  a compact allele -> score table built from the released
  ``UQ_BioSig_model_Final.tsv.gz`` with ``annotatemissense build-score-index``.

Both are keyed by GRCh38 genomic allele. The annotation database may contain
more than one row per allele (one per protein-level annotation join); the
AnnotateMissense score is identical across those rows, and lookups return one
record per allele together with the number of underlying rows.
"""

from __future__ import annotations

import os
import re
import threading
import uuid
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import duckdb

from . import CONFIDENCE_FIELD, PRED_FIELD, SCORE_FIELD

ANNOTATION_TABLE = "variants"
SCORE_TABLE = "scores"
KEY_COLUMNS = ("Chr", "Start", "Ref", "Alt")
MODEL_COLUMNS = ("UQ_score", "UQ_sd", "UQ_pred")

Variant = Tuple[str, int, str, str]

_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}
_ALLELE_RE = re.compile(r"^[ACGTN]+$")

# Columns returned for a single-variant lookup.
DETAIL_COLUMNS: Tuple[str, ...] = (
    "avsnp150", "hg19_chr", "hg19_pos", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene",
    "AAChange.refGene", "HGVS.p", "genename", "protein_variant", "Ensembl_geneid",
    "Ensembl_transcriptid", "GENCODE.v43.ENSG", "GENCODE.v43.ENST", "cytoBand",
    "Interpro_domain", "gnomad_constraint",
    "CLNSIG", "CLNDN", "CLNREVSTAT", "CLNALLELEID", "combined_score",
    "esm1v_t33_650M_UR90S_1", "esm1v_t33_650M_UR90S_2", "esm1v_t33_650M_UR90S_3",
    "esm1v_t33_650M_UR90S_4", "esm1v_t33_650M_UR90S_5", "am_pathogenicity", "am_class",
    "CADD_phred", "REVEL_score", "MetaRNN_score", "MetaSVM_score", "MetaLR_score",
    "BayesDel_addAF_score", "BayesDel_noAF_score", "ClinPred_score", "MVP_score",
    "MutPred_score", "MPC_score", "PrimateAI_score", "DEOGEN2_score", "M-CAP_score",
    "VEST4_score", "DANN_score", "LINSIGHT", "FATHMM_score", "PROVEAN_score", "SIFT_score",
    "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score",
    "LRT_score", "MutationTaster_score", "MutationAssessor_score", "GERP++_RS",
    "phyloP100way_vertebrate", "phastCons100way_vertebrate", "SiPhy_29way_logOdds",
    "integrated_fitCons_score", "gnomad30_AF", "gnomad30_AF_afr", "gnomad30_AF_amr",
    "gnomad30_AF_eas", "gnomad30_AF_nfe", "gnomad30_AF_sas", "gnomad30_AF_fin",
    "AF_popmax", "ExAC_ALL",
)

# Columns returned for batch / VCF / gene lookups.
BATCH_COLUMNS: Tuple[str, ...] = (
    "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "CLNSIG", "am_pathogenicity",
    "am_class", "CADD_phred", "REVEL_score", "MetaRNN_score", "gnomad30_AF",
)

TRANSCRIPT_COLUMNS: Tuple[str, ...] = (
    "genename", "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid", "aaref",
    "aaalt", "aapos", "uniprot_id", "transcript_id", "protein_variant", "GENCODE.v43.ENSG",
    "GENCODE.v43.ENST", "GENCODE.v43.ENSP", "HGVS.p", "AAChange.refGene", "UQ_score",
)

_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


# --------------------------------------------------------------------------- #
# Input normalisation
# --------------------------------------------------------------------------- #
def normalize_chr(chrom: Any) -> str:
    """Return the database chromosome label (``chr1`` ... ``chr22``, ``chrX``, ``chrY``, ``chrM``)."""
    value = str(chrom if chrom is not None else "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "MT":
        value = "M"
    if value.isdigit():
        value = str(int(value))
    if value not in _VALID_CHROMS:
        raise ValueError(f"Unrecognised chromosome: {chrom!r} (expected 1-22, X, Y, M with or without 'chr')")
    return f"chr{value}"


def normalize_variant(chrom: Any, pos: Any, ref: Any, alt: Any) -> Variant:
    """Validate and normalise one GRCh38 allele. Raises ``ValueError`` on bad input."""
    norm_chr = normalize_chr(chrom)
    try:
        norm_pos = int(str(pos).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Position must be an integer: {pos!r}") from exc
    if norm_pos < 1:
        raise ValueError("Position must be 1-based and positive")
    norm_ref = str(ref or "").strip().upper()
    norm_alt = str(alt or "").strip().upper()
    if not norm_ref or not norm_alt:
        raise ValueError("Reference and alternate alleles are required")
    if not _ALLELE_RE.match(norm_ref) or not _ALLELE_RE.match(norm_alt):
        raise ValueError(f"Alleles must contain only A, C, G, T or N: {ref!r}>{alt!r}")
    return norm_chr, norm_pos, norm_ref, norm_alt


def trim_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Reduce REF/ALT to a minimal representation (e.g. ``AT>GT`` -> ``A>G``)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and ref != alt


def parse_variant_string(text: str) -> Variant:
    """Parse ``22:15528192:A:C``, ``chr22-15528192-A-C``, ``22 15528192 A C`` or ``22:15528192 A>C``."""
    parts = [p for p in re.split(r"[\s:,\-_>/]+", text.strip()) if p]
    if len(parts) != 4:
        raise ValueError(f"Cannot parse variant {text!r}; use CHR:POS:REF:ALT")
    return normalize_variant(*parts)


def column_group(name: str) -> str:
    """Assign a released-database column to its annotation source (used by /schema)."""
    if name in KEY_COLUMNS or name == "End":
        return "Genomic allele key (ANNOVAR input, GRCh38)"
    if name in MODEL_COLUMNS:
        return "AnnotateMissense prediction"
    if name.startswith("dbnsfp42c_"):
        return "dbNSFP v4.2c (ANNOVAR filter-based)"
    if name.startswith(("esm1v_", "GENCODE.v43.")) or name in {"HGVS.p", "combined_score"}:
        return "ESM1v variant-effect predictions"
    if name in {"am_pathogenicity", "am_class", "uniprot_id", "transcript_id", "protein_variant", "AlphaMissense"}:
        return "AlphaMissense"
    if name.startswith("CLN"):
        return "ClinVar (ANNOVAR filter-based)"
    if re.match(r"^(Func|Gene|GeneDetail|ExonicFunc|AAChange)\.", name):
        return "ANNOVAR gene-based (" + name.split(".", 1)[1] + ")"
    if name in {"cytoBand", "genomicSuperDups", "rmsk", "simpleRepeat", "phastCons100way"}:
        return "ANNOVAR region-based"
    if name.endswith(".sites.2015_08"):
        return "1000 Genomes (2015-08)"
    if name.startswith("ExAC_"):
        return "ExAC"
    if name.startswith("gnomad211_"):
        return "gnomAD v2.1.1 exome"
    if name.startswith("gnomad30_") or name == "AF_ami":
        return "gnomAD v3.0 genome"
    if name.startswith(("AF", "non_topmed_", "non_neuro_", "non_cancer_", "controls_")):
        return "gnomAD exome (ANNOVAR gnomad_exome)"
    if name in {"avsnp150"}:
        return "dbSNP 150"
    if name in {"cosmic70"}:
        return "COSMIC v70"
    if name.startswith("dbscSNV_"):
        return "dbscSNV splicing"
    if name == "gnomad_constraint":
        return "gnomAD gene constraint"
    if name in {"aaref", "aaalt", "aapos", "genename", "cds_strand", "refcodon", "codonpos",
                "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid", "hg19_chr",
                "hg19_pos", "hg18_chr", "hg18_pos", "genome"}:
        return "dbNSFP v5.1 transcript/protein context (challenge input)"
    return "dbNSFP v4.2a (ANNOVAR filter-based)"


def _qi(identifier: str) -> str:
    """Quote a SQL identifier."""
    return '"' + identifier.replace('"', '""') + '"'


def _ql(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def _split(value: Any, sep: str = ";") -> List[str]:
    if value is None:
        return []
    out: List[str] = []
    for item in str(value).split(sep):
        item = item.strip()
        if item and item != "." and item not in out:
            out.append(item)
    return out


def _hgvs_p_residues(hgvs_p: Optional[str]) -> Optional[Tuple[str, str]]:
    if not hgvs_p:
        return None
    match = re.search(r"p\.\(?([A-Z][a-z]{2})\d+([A-Z][a-z]{2})", hgvs_p)
    if not match:
        return None
    ref, alt = _THREE_TO_ONE.get(match.group(1)), _THREE_TO_ONE.get(match.group(2))
    return (ref, alt) if ref and alt else None


# --------------------------------------------------------------------------- #
# Database access
# --------------------------------------------------------------------------- #
class AnnotateMissenseDB:
    """Thread-safe, read-only access to the annotation database and/or score index."""

    def __init__(self, annotations_path: Optional[str] = None, scores_path: Optional[str] = None) -> None:
        if not annotations_path and not scores_path:
            raise ValueError("Provide an annotation database (variants.duckdb) and/or a score index")
        self.annotations_path = annotations_path or None
        self.scores_path = scores_path or None
        self._conn = duckdb.connect(":memory:")
        self._local = threading.local()
        self.annotation_columns: List[str] = []
        self.annotation_types: Dict[str, str] = {}

        if self.annotations_path:
            self._attach(self.annotations_path, "ann")
            described = self._conn.execute(
                "SELECT column_name, data_type FROM information_schema.columns "
                "WHERE table_catalog = 'ann' AND table_name = ? ORDER BY ordinal_position",
                [ANNOTATION_TABLE],
            ).fetchall()
            if not described:
                raise ValueError(f"{self.annotations_path} has no '{ANNOTATION_TABLE}' table")
            self.annotation_columns = [name for name, _ in described]
            self.annotation_types = dict(described)
            missing = [c for c in KEY_COLUMNS if c not in self.annotation_types]
            if missing:
                raise ValueError(f"'{ANNOTATION_TABLE}' is missing key columns: {missing}")

        if self.scores_path:
            self._attach(self.scores_path, "sc")
            tables = self._conn.execute(
                "SELECT table_name FROM information_schema.tables WHERE table_catalog = 'sc'"
            ).fetchall()
            if (SCORE_TABLE,) not in tables:
                raise ValueError(f"{self.scores_path} has no '{SCORE_TABLE}' table")

        self._column_set = set(self.annotation_columns)

    # -- connection helpers ------------------------------------------------- #
    def _attach(self, path: str, alias: str) -> None:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"DuckDB file not found: {path}")
        self._conn.execute(f"ATTACH {_ql(os.path.abspath(path))} AS {alias} (READ_ONLY)")

    def _cursor(self) -> duckdb.DuckDBPyConnection:
        cur = getattr(self._local, "cursor", None)
        if cur is None:
            cur = self._conn.cursor()
            self._local.cursor = cur
        return cur

    def query(self, sql: str, params: Optional[Sequence[Any]] = None) -> List[Dict[str, Any]]:
        cur = self._cursor().execute(sql, list(params or []))
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]

    @property
    def has_annotations(self) -> bool:
        return bool(self.annotations_path)

    @property
    def has_scores(self) -> bool:
        return bool(self.scores_path)

    def describe(self) -> Dict[str, Any]:
        return {
            "annotation_database": os.path.basename(self.annotations_path) if self.annotations_path else None,
            "score_index": os.path.basename(self.scores_path) if self.scores_path else None,
            "annotation_columns": len(self.annotation_columns),
        }

    # -- core lookup -------------------------------------------------------- #
    def lookup(self, variants: Sequence[Variant], columns: Iterable[str] = BATCH_COLUMNS,
               _query_table: Optional[str] = None) -> List[Dict[str, Any]]:
        """Look up alleles; returns one dict per input variant, in input order.

        Every record contains ``query_*`` keys, ``found``, the canonical
        AnnotateMissense fields and the requested annotation columns (``None``
        when unavailable). ``annotation_rows`` reports how many database rows
        share the allele.
        """
        if not variants:
            return []
        wanted = [c for c in dict.fromkeys(columns) if c not in KEY_COLUMNS]
        present = [c for c in wanted if c in self._column_set]
        absent = [c for c in wanted if c not in self._column_set]

        params: List[Any] = []
        if _query_table is None:
            for idx, (chrom, pos, ref, alt) in enumerate(variants):
                params.extend([idx, chrom, int(pos), ref, alt])
            q_cte = f"q(idx, Chr, Start, Ref, Alt) AS (VALUES {', '.join(['(?, ?, ?, ?, ?)'] * len(variants))})"
        else:
            q_cte = f"q AS (SELECT idx, Chr, Start, Ref, Alt FROM {_query_table})"

        selects = ["q.idx", "q.Chr AS query_chr", "q.Start AS query_pos", "q.Ref AS query_ref", "q.Alt AS query_alt"]
        joins: List[str] = []
        score_terms, conf_terms, pred_terms, found_terms = [], [], [], []

        if self.has_scores:
            joins.append(
                "LEFT JOIN (SELECT s.chrom, s.pos, s.ref, s.alt, s.score, s.confidence, s.pred "
                f"FROM sc.{SCORE_TABLE} s SEMI JOIN q ON s.chrom = q.Chr AND s.pos = q.Start "
                "AND s.ref = q.Ref AND s.alt = q.Alt) s "
                "ON s.chrom = q.Chr AND s.pos = q.Start AND s.ref = q.Ref AND s.alt = q.Alt"
            )
            score_terms.append("s.score")
            conf_terms.append("s.confidence")
            pred_terms.append("s.pred")
            found_terms.append("s.pos IS NOT NULL")
            selects.append("s.pos IS NOT NULL AS in_score_index")

        if self.has_annotations:
            inner = list(dict.fromkeys(list(KEY_COLUMNS) + present + [c for c in MODEL_COLUMNS if c in self._column_set]))
            inner_sql = ", ".join(f"v.{_qi(c)}" for c in inner)
            key = "v.Chr, v.Start, v.Ref, v.Alt"
            joins.append(
                f"LEFT JOIN (SELECT {inner_sql}, COUNT(*) OVER (PARTITION BY {key}) AS _annotation_rows, "
                f"ROW_NUMBER() OVER (PARTITION BY {key} ORDER BY v.rowid) AS _rn "
                f"FROM ann.{ANNOTATION_TABLE} v SEMI JOIN q ON v.Chr = q.Chr AND v.Start = q.Start "
                "AND v.Ref = q.Ref AND v.Alt = q.Alt) a "
                "ON a.Chr = q.Chr AND a.Start = q.Start AND a.Ref = q.Ref AND a.Alt = q.Alt AND a._rn = 1"
            )
            selects += [f"a.{_qi(c)}" for c in present]
            selects.append("a._annotation_rows AS annotation_rows")
            if "UQ_score" in self._column_set:
                score_terms.append("a.UQ_score")
            if "UQ_sd" in self._column_set:
                conf_terms.append("a.UQ_sd")
            if "UQ_pred" in self._column_set:
                pred_terms.append("a.UQ_pred")
            found_terms.append("a._annotation_rows IS NOT NULL")
        else:
            selects.append("NULL AS annotation_rows")

        selects += [f"NULL AS {_qi(c)}" for c in absent]

        def coalesce(terms: List[str]) -> str:
            return f"COALESCE({', '.join(terms)})" if terms else "NULL"

        selects += [
            f"{coalesce(score_terms)} AS {SCORE_FIELD}",
            f"{coalesce(conf_terms)} AS {CONFIDENCE_FIELD}",
            f"{coalesce(pred_terms)} AS {PRED_FIELD}",
            f"({' OR '.join(found_terms)}) AS found",
        ]

        sql = (
            f"WITH {q_cte} "
            f"SELECT {', '.join(selects)} FROM q {' '.join(joins)} ORDER BY q.idx"
        )
        rows = self.query(sql, params)
        for row in rows:
            row["found"] = bool(row["found"])
            if row.get("annotation_rows") is not None:
                row["annotation_rows"] = int(row["annotation_rows"])
        return rows

    def lookup_chunked(self, variants: Sequence[Variant], columns: Iterable[str] = BATCH_COLUMNS,
                       chunk_size: int = 1000) -> List[Dict[str, Any]]:
        """Lookup for any number of alleles: small inputs inline, large inputs via one temp-table join
        (a single scan of the database instead of one scan per chunk)."""
        columns = tuple(columns)
        if len(variants) <= chunk_size:
            return self.lookup(variants, columns)
        cur = self._cursor()
        name = f"_am_query_{uuid.uuid4().hex}"
        cur.execute(f"CREATE TEMP TABLE {name} (idx BIGINT, Chr VARCHAR, Start BIGINT, Ref VARCHAR, Alt VARCHAR)")
        try:
            for start in range(0, len(variants), chunk_size):
                chunk = variants[start:start + chunk_size]
                params: List[Any] = []
                for offset, (chrom, pos, ref, alt) in enumerate(chunk):
                    params.extend([start + offset, chrom, int(pos), ref, alt])
                cur.execute(f"INSERT INTO {name} VALUES " + ", ".join(["(?, ?, ?, ?, ?)"] * len(chunk)), params)
            return self.lookup(variants, columns, _query_table=name)
        finally:
            cur.execute(f"DROP TABLE IF EXISTS {name}")

    def fetch_variant(self, variant: Variant) -> Optional[Dict[str, Any]]:
        row = self.lookup([variant], DETAIL_COLUMNS)[0]
        return row if row["found"] else None

    def fetch_variant_full(self, variant: Variant) -> Optional[Dict[str, Any]]:
        """All 337 columns of the first database row for the allele (annotation DB only)."""
        if not self.has_annotations:
            return None
        rows = self.query(
            f"SELECT * FROM ann.{ANNOTATION_TABLE} WHERE Chr = ? AND Start = ? AND Ref = ? AND Alt = ? "
            "ORDER BY rowid LIMIT 1",
            list(variant),
        )
        return rows[0] if rows else None

    # -- transcript provenance ---------------------------------------------- #
    def transcript_context(self, variant: Variant) -> Optional[Dict[str, Any]]:
        """Explain which transcript/protein each annotation source used for this allele."""
        if not self.has_annotations:
            return None
        cols = [c for c in TRANSCRIPT_COLUMNS if c in self._column_set]
        rows = self.query(
            f"SELECT {', '.join(_qi(c) for c in cols)} FROM ann.{ANNOTATION_TABLE} "
            "WHERE Chr = ? AND Start = ? AND Ref = ? AND Alt = ? ORDER BY rowid",
            list(variant),
        )
        if not rows:
            return None

        def union(col: str, sep: str = ";") -> List[str]:
            out: List[str] = []
            for r in rows:
                for item in _split(r.get(col), sep):
                    if item not in out:
                        out.append(item)
            return out

        first = rows[0]
        aaref, aaalt = first.get("aaref"), first.get("aaalt")
        am_variant = first.get("protein_variant")
        am_match = None
        if am_variant and aaref and aaalt:
            am_match = am_variant[0] == aaref and am_variant[-1] == aaalt

        esm_records: List[Dict[str, Any]] = []
        seen = set()
        for r in rows:
            key = (r.get("GENCODE.v43.ENST"), r.get("HGVS.p"))
            if key == (None, None) or key in seen:
                continue
            seen.add(key)
            residues = _hgvs_p_residues(r.get("HGVS.p"))
            esm_records.append({
                "transcript": r.get("GENCODE.v43.ENST"),
                "protein": r.get("GENCODE.v43.ENSP"),
                "hgvs_p": r.get("HGVS.p"),
                "matches_dbnsfp_amino_acids": (residues == (aaref, aaalt)) if residues and aaref else None,
            })

        scores = {r.get("UQ_score") for r in rows}
        return {
            "annotation_rows": len(rows),
            "distinct_annotatemissense_scores": len(scores),
            "annotatemissense_score_is_per_allele": len(scores) <= 1,
            "dbnsfp": {
                "genes": union("genename"),
                "ensembl_gene": union("Ensembl_geneid"),
                "transcripts": union("Ensembl_transcriptid"),
                "proteins": union("Ensembl_proteinid"),
                "aa_positions": union("aapos"),
                "ref_aa": aaref,
                "alt_aa": aaalt,
            },
            "alphamissense": {
                "transcript": first.get("transcript_id"),
                "uniprot": first.get("uniprot_id"),
                "protein_variant": am_variant,
                "matches_dbnsfp_amino_acids": am_match,
            },
            "esm1v": esm_records,
            "annovar_refgene_aachange": union("AAChange.refGene", ","),
        }

    # -- gene / stats / schema ---------------------------------------------- #
    def gene_variants(self, gene: str, limit: int = 100) -> List[Dict[str, Any]]:
        if not self.has_annotations:
            return []
        gene_l = gene.strip().lower()
        keys = self.query(
            f"SELECT Chr, Start, Ref, Alt FROM ann.{ANNOTATION_TABLE} "
            "WHERE list_contains(string_split(lower(\"Gene.refGene\"), ';'), ?) "
            "   OR list_contains(string_split(lower(genename), ';'), ?) "
            "GROUP BY Chr, Start, Ref, Alt ORDER BY Chr, Start, Ref, Alt LIMIT ?",
            [gene_l, gene_l, int(limit)],
        )
        variants = [(k["Chr"], int(k["Start"]), k["Ref"], k["Alt"]) for k in keys]
        return self.lookup(variants, BATCH_COLUMNS) if variants else []

    def stats(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"resources": self.describe()}
        if self.has_scores:
            row = self.query(
                f"SELECT COUNT(*) AS alleles, COUNT(score) AS scored, "
                f"SUM(CASE WHEN pred = 'D' THEN 1 ELSE 0 END) AS damaging, "
                f"SUM(CASE WHEN pred = 'T' THEN 1 ELSE 0 END) AS tolerated FROM sc.{SCORE_TABLE}"
            )[0]
            per_chr = self.query(f"SELECT chrom, COUNT(*) AS n FROM sc.{SCORE_TABLE} GROUP BY 1 ORDER BY 1")
            out["score_index"] = {
                "alleles": int(row["alleles"]),
                "scored_alleles": int(row["scored"]),
                "damaging": int(row["damaging"] or 0),
                "tolerated": int(row["tolerated"] or 0),
                "chromosomes": {r["chrom"]: int(r["n"]) for r in per_chr},
            }
        if self.has_annotations:
            row = self.query(
                f"SELECT COUNT(*) AS n_rows, COUNT(CLNSIG) AS with_clinvar, COUNT(UQ_score) AS with_score, "
                f"SUM(CASE WHEN UQ_pred = 'D' THEN 1 ELSE 0 END) AS damaging, "
                f"SUM(CASE WHEN UQ_pred = 'T' THEN 1 ELSE 0 END) AS tolerated FROM ann.{ANNOTATION_TABLE}"
            )[0]
            per_chr = self.query(f"SELECT Chr, COUNT(*) AS n FROM ann.{ANNOTATION_TABLE} GROUP BY 1 ORDER BY 1")
            out["annotation_database"] = {
                "rows": int(row["n_rows"]),
                "columns": len(self.annotation_columns),
                "rows_with_clinvar": int(row["with_clinvar"]),
                "rows_with_score": int(row["with_score"]),
                "damaging_rows": int(row["damaging"] or 0),
                "tolerated_rows": int(row["tolerated"] or 0),
                "chromosomes": {r["Chr"]: int(r["n"]) for r in per_chr},
            }
        return out

    def schema(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        if self.has_annotations:
            out[ANNOTATION_TABLE] = [
                {"position": i + 1, "name": c, "type": self.annotation_types[c], "source": column_group(c)}
                for i, c in enumerate(self.annotation_columns)
            ]
        if self.has_scores:
            rows = self.query(
                "SELECT column_name, data_type FROM information_schema.columns "
                "WHERE table_catalog = 'sc' AND table_name = ? ORDER BY ordinal_position",
                [SCORE_TABLE],
            )
            out[SCORE_TABLE] = [
                {"position": i + 1, "name": r["column_name"], "type": r["data_type"]} for i, r in enumerate(rows)
            ]
        return out

    def audit_alleles(self) -> Dict[str, Any]:
        """Row/allele multiplicity audit of the annotation database (full scan)."""
        if not self.has_annotations:
            raise ValueError("Allele audit requires the annotation database")
        row = self.query(
            f"""
            WITH g AS (
                SELECT Chr, Start, Ref, Alt, COUNT(*) AS n,
                       COUNT(DISTINCT UQ_score) AS n_scores
                FROM ann.{ANNOTATION_TABLE} GROUP BY ALL
            )
            SELECT SUM(n) AS rows, COUNT(*) AS alleles,
                   SUM(CASE WHEN n > 1 THEN 1 ELSE 0 END) AS alleles_with_multiple_rows,
                   MAX(n) AS max_rows_per_allele,
                   SUM(CASE WHEN n_scores > 1 THEN 1 ELSE 0 END) AS alleles_with_discordant_scores
            FROM g
            """
        )[0]
        return {k: int(v or 0) for k, v in row.items()}


def build_score_index(tsv_path: str, out_path: str, threads: Optional[int] = None) -> Dict[str, Any]:
    """Convert ``UQ_BioSig_model_Final.tsv.gz`` into a compact, sorted DuckDB score index.

    * ``pred = 'U'`` rows ("No prediction available") are kept with NULL score.
    * Exact duplicate allele rows are collapsed to one row per allele.
    * ``sd`` is stored as ``confidence``: it equals the probability of the
      predicted class, max(score, 1 - score), not a standard deviation.
    """
    if os.path.exists(out_path):
        raise FileExistsError(f"{out_path} already exists")
    con = duckdb.connect(out_path)
    try:
        if threads:
            con.execute(f"SET threads = {int(threads)}")
        con.execute(
            f"""
            CREATE TABLE {SCORE_TABLE} AS
            WITH raw AS (
                SELECT
                    'chr' || CASE WHEN upper(CAST("#chr" AS VARCHAR)) IN ('MT', 'M') THEN 'M'
                                  ELSE upper(CAST("#chr" AS VARCHAR)) END AS chrom,
                    CAST("pos(1-based)" AS INTEGER) AS pos,
                    upper(ref) AS ref, upper(alt) AS alt,
                    CASE WHEN pred IN ('D', 'T') THEN CAST(score AS FLOAT) END AS score,
                    CASE WHEN pred IN ('D', 'T') THEN CAST(sd AS FLOAT) END AS confidence,
                    CASE WHEN pred IN ('D', 'T') THEN pred END AS pred
                FROM read_csv({_ql(os.path.abspath(tsv_path))}, delim = '\t', header = true,
                              columns = {{'#chr': 'VARCHAR', 'pos(1-based)': 'BIGINT', 'ref': 'VARCHAR',
                                         'alt': 'VARCHAR', 'score': 'DOUBLE', 'sd': 'DOUBLE',
                                         'pred': 'VARCHAR', 'comments': 'VARCHAR'}})
            )
            SELECT chrom, pos, ref, alt,
                   any_value(score) AS score, any_value(confidence) AS confidence, any_value(pred) AS pred,
                   COUNT(*)::USMALLINT AS source_rows
            FROM raw GROUP BY chrom, pos, ref, alt
            ORDER BY chrom, pos, ref, alt
            """
        )
        stats = con.execute(
            f"SELECT SUM(source_rows), COUNT(*), COUNT(score), "
            f"SUM(CASE WHEN source_rows > 1 THEN 1 ELSE 0 END) FROM {SCORE_TABLE}"
        ).fetchone()
        con.execute("CREATE TABLE meta (key VARCHAR, value VARCHAR)")
        meta = {
            "source_file": os.path.basename(tsv_path),
            "source_rows": str(int(stats[0] or 0)),
            "alleles": str(int(stats[1] or 0)),
            "scored_alleles": str(int(stats[2] or 0)),
            "duplicated_alleles_collapsed": str(int(stats[3] or 0)),
            "genome_build": "GRCh38",
        }
        con.executemany("INSERT INTO meta VALUES (?, ?)", list(meta.items()))
    finally:
        con.close()
    meta["output_bytes"] = str(os.path.getsize(out_path))
    return meta
