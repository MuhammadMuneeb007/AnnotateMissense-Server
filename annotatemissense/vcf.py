"""Annotate user-supplied VCFs or variant lists by lookup against the release.

AnnotateMissense does **not** run the model on new alleles. Each ALT allele of
each input record is normalised (chromosome naming, upper-case, minimal
REF/ALT trimming) and looked up in the precomputed GRCh38 resource. Every ALT
receives one of these statuses:

``found``           score available
``no_score``        allele present in the release but without a prediction
``not_in_release``  a valid SNV that is not among the precomputed missense alleles
``not_snv``         indel, MNV, symbolic or missing ALT (not scored by design)
``invalid``         unparseable chromosome/position/allele
"""

from __future__ import annotations

import datetime as _dt
import gzip
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import CONFIDENCE_FIELD, PRED_FIELD, SCORE_FIELD, __version__
from .db import BATCH_COLUMNS, AnnotateMissenseDB, Variant, is_snv, normalize_chr, trim_alleles

INFO_SCORE = "AnnotateMissense_score"
INFO_PRED = "AnnotateMissense_pred"
INFO_STATUS = "AnnotateMissense_status"

INFO_HEADER_LINES = [
    f'##INFO=<ID={INFO_SCORE},Number=A,Type=Float,Description="AnnotateMissense pathogenicity score '
    '(0-1; >=0.5 predicted damaging). Precomputed GRCh38 lookup; . when unavailable">',
    f'##INFO=<ID={INFO_PRED},Number=A,Type=String,Description="AnnotateMissense prediction: '
    'D = damaging (score >= 0.5), T = tolerated">',
    f'##INFO=<ID={INFO_STATUS},Number=A,Type=String,Description="AnnotateMissense lookup status: '
    'found | no_score | not_in_release | not_snv | invalid">',
]

TSV_COLUMNS = [
    "line", "chrom", "pos", "id", "ref", "alt", "query", "status", SCORE_FIELD, CONFIDENCE_FIELD,
    PRED_FIELD, "gene", "exonic_function", "aa_change_refgene", "clinvar_clnsig", "alphamissense",
    "cadd_phred", "revel", "metarnn", "gnomad30_af", "annotation_rows",
]

_GRCH37_HINT = re.compile(r"grch37|hg19|\bb37\b|hs37|human_g1k_v37", re.IGNORECASE)
_HEADER_WORDS = {"chr", "chrom", "chromosome", "#chr", "#chrom"}


def decode_upload(data: bytes) -> str:
    """Decode plain or gzip/bgzip-compressed text."""
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


@dataclass
class Entry:
    """One ALT allele of one input record."""

    line_no: int
    chrom: str
    pos: str
    id: str
    ref: str
    alt: str
    status: str = "pending"
    variant: Optional[Variant] = None
    message: Optional[str] = None
    result: Dict[str, Any] = field(default_factory=dict)

    @property
    def query(self) -> str:
        return ":".join(str(x) for x in self.variant) if self.variant else ""

    def as_dict(self) -> Dict[str, Any]:
        r = self.result
        return {
            "line": self.line_no, "chrom": self.chrom, "pos": self.pos, "id": self.id,
            "ref": self.ref, "alt": self.alt, "query": self.query or None, "status": self.status,
            "message": self.message,
            SCORE_FIELD: r.get(SCORE_FIELD), CONFIDENCE_FIELD: r.get(CONFIDENCE_FIELD),
            PRED_FIELD: r.get(PRED_FIELD), "gene": r.get("Gene.refGene"),
            "exonic_function": r.get("ExonicFunc.refGene"), "aa_change_refgene": r.get("AAChange.refGene"),
            "clinvar_clnsig": r.get("CLNSIG"), "alphamissense": r.get("am_pathogenicity"),
            "cadd_phred": r.get("CADD_phred"), "revel": r.get("REVEL_score"),
            "metarnn": r.get("MetaRNN_score"), "gnomad30_af": r.get("gnomad30_AF"),
            "annotation_rows": r.get("annotation_rows"),
        }


@dataclass
class Record:
    line_no: int
    fields: Optional[List[str]]  # original VCF columns (None for simple lists)
    entries: List[Entry]


@dataclass
class AnnotationResult:
    is_vcf: bool
    header: List[str]
    records: List[Record]
    warnings: List[str]

    @property
    def entries(self) -> List[Entry]:
        return [e for rec in self.records for e in rec.entries]

    def summary(self) -> Dict[str, int]:
        counts = {"records": len(self.records), "alt_alleles": 0, "found": 0, "no_score": 0,
                  "not_in_release": 0, "not_snv": 0, "invalid": 0, "damaging": 0, "tolerated": 0}
        for e in self.entries:
            counts["alt_alleles"] += 1
            counts[e.status] = counts.get(e.status, 0) + 1
            pred = e.result.get(PRED_FIELD)
            if pred == "D":
                counts["damaging"] += 1
            elif pred == "T":
                counts["tolerated"] += 1
        return counts


def _split_fields(line: str) -> List[str]:
    """Split a variant-list line: tab, whitespace, comma, or CHR:POS:REF:ALT notation."""
    if "\t" in line:
        return [p.strip() for p in line.split("\t")]
    parts = line.split()
    if len(parts) >= 4:
        return parts
    comma = [p.strip() for p in line.split(",")]
    if len(comma) >= 4:
        return comma
    # 22:36305985:G:A, chr22-36305985-G-A, 22:36305985 G>A (multi-allelic ALT keeps its commas)
    return [p for p in re.split(r"[\s:\-_>/]+", line.strip()) if p]


def _prepare(entry: Entry) -> None:
    alt = entry.alt.strip()
    if not alt or alt in {".", "*"} or alt.startswith("<") or "[" in alt or "]" in alt:
        entry.status = "not_snv"
        entry.message = "Missing, symbolic or breakend ALT"
        return
    try:
        chrom = normalize_chr(entry.chrom)
        pos = int(entry.pos)
        if pos < 1:
            raise ValueError("Position must be positive")
        ref, alt = entry.ref.strip().upper(), alt.upper()
        if not re.match(r"^[ACGTN]+$", ref) or not re.match(r"^[ACGTN]+$", alt):
            raise ValueError(f"Unsupported allele characters: {entry.ref}>{entry.alt}")
    except ValueError as exc:
        entry.status = "invalid"
        entry.message = str(exc)
        return
    pos, ref, alt = trim_alleles(pos, ref, alt)
    if not is_snv(ref, alt):
        entry.status = "not_snv"
        entry.message = "Only single-nucleotide missense alleles are precomputed"
        return
    entry.variant = (chrom, pos, ref, alt)


def parse(text: str, max_alleles: Optional[int] = None) -> AnnotationResult:
    """Parse a VCF or a simple ``chrom pos ref alt`` list (tab, comma or space separated)."""
    header: List[str] = []
    records: List[Record] = []
    is_vcf = False
    n_alleles = 0
    first_data = True

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            header.append(line)
            if line.startswith("##fileformat=VCF") or line.upper().startswith("#CHROM"):
                is_vcf = True
            continue

        if is_vcf:
            fields = line.split("\t")
            if len(fields) < 5:
                fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"Line {line_no}: a VCF record needs at least CHROM POS ID REF ALT")
            chrom, pos, vid, ref, alt_field = fields[:5]
            alts = alt_field.split(",")
            keep_fields: Optional[List[str]] = fields
        else:
            parts = _split_fields(line)
            if first_data and parts and parts[0].lower() in _HEADER_WORDS:
                first_data = False
                continue
            first_data = False
            if len(parts) >= 8 and parts[1].isdigit() and not parts[2].isdigit():
                chrom, pos, vid, ref, alt_field = parts[:5]  # headerless VCF-like row
            elif len(parts) >= 4:
                chrom, pos, ref, alt_field = parts[:4]
                vid = "."
            else:
                raise ValueError(f"Line {line_no}: expected 'chrom pos ref alt'")
            alts = alt_field.split(",")
            keep_fields = None

        entries = [Entry(line_no, chrom, pos, vid, ref, alt) for alt in alts]
        for entry in entries:
            _prepare(entry)
        n_alleles += len(entries)
        if max_alleles is not None and n_alleles > max_alleles:
            raise ValueError(
                f"Input exceeds the limit of {max_alleles:,} ALT alleles for this service; "
                "use the offline 'annotatemissense annotate-vcf' command for larger files"
            )
        records.append(Record(line_no, keep_fields, entries))

    warnings: List[str] = []
    if any(_GRCH37_HINT.search(h) for h in header if h.startswith(("##reference", "##contig", "##assembly"))):
        warnings.append(
            "The VCF header refers to GRCh37/hg19. AnnotateMissense uses GRCh38/hg38 coordinates; "
            "lift the VCF over to GRCh38 first or matches will be wrong or missing."
        )
    return AnnotationResult(is_vcf=is_vcf, header=header, records=records, warnings=warnings)


def annotate(text: str, db: AnnotateMissenseDB, max_alleles: Optional[int] = None,
             columns: Sequence[str] = BATCH_COLUMNS) -> AnnotationResult:
    result = parse(text, max_alleles=max_alleles)
    unique: Dict[Variant, None] = {}
    for e in result.entries:
        if e.variant:
            unique.setdefault(e.variant, None)
    keys = list(unique)
    rows = db.lookup_chunked(keys, columns) if keys else []
    by_variant = {keys[int(r["idx"])]: r for r in rows}
    for e in result.entries:
        if not e.variant:
            continue
        row = by_variant.get(e.variant)
        if row is None or not row["found"]:
            e.status = "not_in_release"
        else:
            e.result = row
            e.status = "found" if row.get(SCORE_FIELD) is not None else "no_score"
    return result


def _fmt(value: Any) -> str:
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def to_vcf(result: AnnotationResult) -> str:
    """Return a VCF with AnnotateMissense INFO fields (Number=A, one value per ALT)."""
    meta = [h for h in result.header if h.startswith("##")]
    column_header = next((h for h in result.header if h.upper().startswith("#CHROM")), None)
    meta = [h for h in meta if not any(h.startswith(f"##INFO=<ID={k},") for k in (INFO_SCORE, INFO_PRED, INFO_STATUS))]
    if not meta or not meta[0].startswith("##fileformat"):
        meta.insert(0, "##fileformat=VCFv4.2")
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    meta += INFO_HEADER_LINES
    meta.append(f'##AnnotateMissense=<Version="{__version__}",Build="GRCh38",Date="{stamp}",'
                'Source="https://doi.org/10.5281/zenodo.19981867">')
    lines = meta + [column_header or "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]

    for rec in result.records:
        if rec.fields is not None:
            fields = list(rec.fields) + ["."] * max(0, 8 - len(rec.fields))
        else:
            e0 = rec.entries[0]
            fields = [e0.chrom, e0.pos, e0.id, e0.ref, ",".join(e.alt for e in rec.entries), ".", ".", "."]
        info = [kv for kv in fields[7].split(";") if kv and kv != "."
                and kv.split("=", 1)[0] not in (INFO_SCORE, INFO_PRED, INFO_STATUS)]
        info.append(f"{INFO_SCORE}=" + ",".join(_fmt(e.result.get(SCORE_FIELD)) for e in rec.entries))
        info.append(f"{INFO_PRED}=" + ",".join(_fmt(e.result.get(PRED_FIELD)) for e in rec.entries))
        info.append(f"{INFO_STATUS}=" + ",".join(e.status for e in rec.entries))
        fields[7] = ";".join(info)
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def to_tsv(result: AnnotationResult) -> str:
    lines = ["\t".join(TSV_COLUMNS)]
    for e in result.entries:
        d = e.as_dict()
        lines.append("\t".join("" if d.get(c) is None else _fmt(d.get(c)) for c in TSV_COLUMNS))
    return "\n".join(lines) + "\n"


def to_json(result: AnnotationResult, limit: Optional[int] = None) -> Dict[str, Any]:
    entries = result.entries
    return {
        "status": "success",
        "input_type": "vcf" if result.is_vcf else "variant_list",
        "genome_build": "GRCh38",
        "summary": result.summary(),
        "warnings": result.warnings,
        "results": [e.as_dict() for e in (entries if limit is None else entries[:limit])],
        "truncated": limit is not None and len(entries) > limit,
    }


def variants_from_text(text: str, max_alleles: Optional[int] = None) -> Tuple[List[Variant], List[Entry]]:
    """Parse text and return (queryable variants, all entries)."""
    result = parse(text, max_alleles=max_alleles)
    entries = result.entries
    return [e.variant for e in entries if e.variant], entries
