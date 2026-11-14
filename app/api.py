"""REST API: /api/v1/*  (JSON; CORS enabled; no cookies; nothing stored)."""

import json
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from flask import Blueprint, Response, current_app, jsonify, request, send_from_directory

from annotatemissense import (CONFIDENCE_FIELD, DECISION_THRESHOLD, GENOME_BUILD, PRED_FIELD,
                              SCORE_FIELD, ZENODO_DOI)
from annotatemissense.db import normalize_variant, parse_variant_string
from annotatemissense.vcf import annotate, decode_upload, parse, to_json, to_tsv, to_vcf

from . import limiter

api_bp = Blueprint("api", __name__, url_prefix="/api/v1")
_CACHE: Dict[str, Any] = {}
_CACHE_TTL = 3600

PRED_LABELS = {"D": "Damaging", "T": "Tolerated"}


def _limit() -> str:
    return current_app.config["RATE_LIMIT"]


def _db():
    return current_app.config["VARIANT_DB"]


def _error(message: str, code: int = 400):
    return jsonify({"status": "error", "message": message}), code


def _param(payload: Dict[str, Any], *names: str) -> Optional[Any]:
    for name in names:
        for source in (request.form, payload, request.args):
            value = source.get(name) if hasattr(source, "get") else None
            if value not in (None, ""):
                return value
    return None


def _variant_from_request(spec: Optional[str] = None):
    if spec:
        return parse_variant_string(spec)
    payload = request.get_json(silent=True) or {}
    chrom = _param(payload, "chrom", "chr")
    pos = _param(payload, "pos", "position")
    ref = _param(payload, "ref")
    alt = _param(payload, "alt")
    if not all([chrom, pos, ref, alt]):
        raise ValueError("Required fields: chrom, pos, ref, alt (or /api/v1/variant/<chrom>:<pos>:<ref>:<alt>)")
    return normalize_variant(chrom, pos, ref, alt)


def _prediction_block(row: Dict[str, Any]) -> Dict[str, Any]:
    pred = row.get(PRED_FIELD)
    return {
        "score": row.get(SCORE_FIELD),
        "confidence": row.get(CONFIDENCE_FIELD),
        "prediction": pred,
        "label": PRED_LABELS.get(pred, "No prediction"),
        "threshold": DECISION_THRESHOLD,
        "note": "confidence = probability of the predicted class, max(score, 1 - score)",
    }


_INT_TEXT = re.compile(r"^-?\d+\.0$")


def _clean(value: Any) -> Any:
    """Tidy text artefacts of the source annotations ('36702031.0', ANNOVAR '\\x2c' escapes)."""
    if isinstance(value, str):
        if _INT_TEXT.match(value):
            return value[:-2]
        return value.replace("\\x2c", ",").replace("\\x3b", ";").replace("\\x3d", "=")
    return value


def format_single_variant(row: Dict[str, Any], transcripts: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    def g(key: str) -> Any:
        return _clean(row.get(key))
    return {
        "status": "success",
        "genome_build": GENOME_BUILD,
        "variant": {
            "chr": g("query_chr"), "pos": g("query_pos"), "ref": g("query_ref"), "alt": g("query_alt"),
            "rsid": g("avsnp150"), "hgvs_p": g("HGVS.p") or g("protein_variant"),
            "hg19": {"chr": g("hg19_chr"), "pos": g("hg19_pos")},
        },
        "annotatemissense": _prediction_block(row),
        "gene": {
            "name": g("Gene.refGene") or g("genename"), "region": g("Func.refGene"),
            "function": g("ExonicFunc.refGene"), "aa_change": g("AAChange.refGene"),
            "cytoband": g("cytoBand"), "ensembl_gene": g("Ensembl_geneid") or g("GENCODE.v43.ENSG"),
            "ensembl_transcript": g("Ensembl_transcriptid") or g("GENCODE.v43.ENST"),
            "interpro_domain": g("Interpro_domain"), "gnomad_constraint": g("gnomad_constraint"),
        },
        "clinical": {"CLNSIG": g("CLNSIG"), "CLNDN": g("CLNDN"), "CLNREVSTAT": g("CLNREVSTAT"),
                     "CLNALLELEID": g("CLNALLELEID")},
        "esm1v": {
            "models": [g(f"esm1v_t33_650M_UR90S_{i}") for i in range(1, 6)],
            "combined_score": g("combined_score"),
        },
        "pathogenicity_scores": {
            "AlphaMissense": g("am_pathogenicity"), "AlphaMissense_class": g("am_class"),
            "CADD_phred": g("CADD_phred"), "REVEL_score": g("REVEL_score"), "MetaRNN_score": g("MetaRNN_score"),
            "MetaSVM_score": g("MetaSVM_score"), "MetaLR_score": g("MetaLR_score"),
            "BayesDel_addAF_score": g("BayesDel_addAF_score"), "BayesDel_noAF_score": g("BayesDel_noAF_score"),
            "ClinPred_score": g("ClinPred_score"), "MVP_score": g("MVP_score"), "MutPred_score": g("MutPred_score"),
            "MPC_score": g("MPC_score"), "PrimateAI_score": g("PrimateAI_score"), "DEOGEN2_score": g("DEOGEN2_score"),
            "M-CAP_score": g("M-CAP_score"), "VEST4_score": g("VEST4_score"), "DANN_score": g("DANN_score"),
            "LINSIGHT": g("LINSIGHT"), "FATHMM_score": g("FATHMM_score"), "PROVEAN_score": g("PROVEAN_score"),
            "SIFT_score": g("SIFT_score"), "SIFT_pred": g("SIFT_pred"),
            "Polyphen2_HDIV_score": g("Polyphen2_HDIV_score"), "Polyphen2_HDIV_pred": g("Polyphen2_HDIV_pred"),
            "Polyphen2_HVAR_score": g("Polyphen2_HVAR_score"), "LRT_score": g("LRT_score"),
            "MutationTaster_score": g("MutationTaster_score"), "MutationAssessor_score": g("MutationAssessor_score"),
        },
        "conservation": {
            "GERP++_RS": g("GERP++_RS"), "phyloP100way_vertebrate": g("phyloP100way_vertebrate"),
            "phastCons100way_vertebrate": g("phastCons100way_vertebrate"),
            "SiPhy_29way_logOdds": g("SiPhy_29way_logOdds"), "integrated_fitCons_score": g("integrated_fitCons_score"),
        },
        "population": {
            "gnomad30_AF": g("gnomad30_AF"), "AF_popmax": g("AF_popmax"), "gnomad30_AF_afr": g("gnomad30_AF_afr"),
            "gnomad30_AF_amr": g("gnomad30_AF_amr"), "gnomad30_AF_eas": g("gnomad30_AF_eas"),
            "gnomad30_AF_nfe": g("gnomad30_AF_nfe"), "gnomad30_AF_sas": g("gnomad30_AF_sas"),
            "gnomad30_AF_fin": g("gnomad30_AF_fin"), "ExAC_ALL": g("ExAC_ALL"),
        },
        "transcripts": transcripts,
        "provenance": {
            "annotation_rows": g("annotation_rows"),
            "in_score_index": g("in_score_index"),
            "source": f"https://doi.org/{ZENODO_DOI}",
        },
    }


def _batch_row(row: Dict[str, Any]) -> Dict[str, Any]:
    chrom, pos, ref, alt = row.get("query_chr"), row.get("query_pos"), row.get("query_ref"), row.get("query_alt")
    return {
        "chr": chrom, "pos": pos, "ref": ref, "alt": alt, "variant": f"{chrom}:{pos} {ref}>{alt}",
        "found": row.get("found"),
        "status": ("found" if row.get(SCORE_FIELD) is not None else "no_score") if row.get("found") else "not_in_release",
        "gene": row.get("Gene.refGene"), "exonic_func": row.get("ExonicFunc.refGene"),
        "aa_change": row.get("AAChange.refGene"), "CLNSIG": _clean(row.get("CLNSIG")),
        SCORE_FIELD: row.get(SCORE_FIELD), CONFIDENCE_FIELD: row.get(CONFIDENCE_FIELD), PRED_FIELD: row.get(PRED_FIELD),
        "AlphaMissense": row.get("am_pathogenicity"), "AlphaMissense_class": row.get("am_class"),
        "CADD_phred": row.get("CADD_phred"), "REVEL_score": row.get("REVEL_score"),
        "MetaRNN_score": row.get("MetaRNN_score"), "gnomad30_AF": row.get("gnomad30_AF"),
        "annotation_rows": row.get("annotation_rows"),
    }


# --------------------------------------------------------------------------- #
@api_bp.get("/health")
def api_health():
    db = _db()
    cov = current_app.config["COVERAGE"]
    return jsonify({
        "status": "ok", "version": current_app.config["VERSION"], "genome_build": GENOME_BUILD,
        "resources": db.describe(), "coverage": cov,
        "limits": {"batch_variants": current_app.config["MAX_BATCH"],
                   "vcf_alt_alleles": current_app.config["MAX_VCF_ALLELES"],
                   "rate_limit": current_app.config["RATE_LIMIT"]},
    })


@api_bp.route("/variant", methods=["GET", "POST"])
@api_bp.get("/variant/<path:spec>")
@limiter.limit(_limit)
def api_variant(spec: Optional[str] = None):
    try:
        variant = _variant_from_request(spec)
    except ValueError as exc:
        return _error(str(exc), 400)
    db = _db()
    row = db.fetch_variant(variant)
    if row is None:
        return jsonify({
            "status": "error", "message": "Variant not found",
            "detail": "Allele is not among the precomputed GRCh38 missense SNVs served by this instance. "
                      "Absence is not evidence of benignity.",
            "query": dict(zip(("chr", "pos", "ref", "alt"), variant)),
            "coverage": current_app.config["COVERAGE"]["label"],
        }), 404
    return jsonify(format_single_variant(row, db.transcript_context(variant)))


@api_bp.route("/variant/full", methods=["GET", "POST"])
@limiter.limit(_limit)
def api_variant_full():
    try:
        variant = _variant_from_request()
    except ValueError as exc:
        return _error(str(exc), 400)
    row = _db().fetch_variant_full(variant)
    if row is None:
        return _error("Variant not found in the annotation database", 404)
    return jsonify({
        "status": "success", "genome_build": GENOME_BUILD,
        "variant": dict(zip(("chr", "pos", "ref", "alt"), variant)),
        "note": "First database row for the allele; legacy columns UQ_score/UQ_sd/UQ_pred are the "
                "AnnotateMissense score/confidence/prediction.",
        "all_columns": row,
    })


@api_bp.post("/batch")
@limiter.limit(_limit)
def api_batch():
    texts: List[str] = []
    upload = request.files.get("variants_file") or request.files.get("file")
    if upload and upload.filename:
        texts.append(decode_upload(upload.read()))
    payload = request.get_json(silent=True) or {}
    if request.form.get("variants_text", "").strip():
        texts.append(request.form["variants_text"])
    if isinstance(payload.get("variants"), list):
        texts.append("\n".join(
            " ".join(str(v.get(k, "")) for k in ("chr", "pos", "ref", "alt")) if isinstance(v, dict) else str(v)
            for v in payload["variants"]
        ))
    elif isinstance(payload.get("variants"), str):
        texts.append(payload["variants"])
    if not texts:
        return _error("Provide variants_file, variants_text, or JSON {\"variants\": [...]}", 400)

    max_batch = current_app.config["MAX_BATCH"]
    try:
        entries = [e for t in texts for e in parse(t, max_alleles=max_batch).entries]
    except ValueError as exc:
        return _error(str(exc), 400)
    if len(entries) > max_batch:
        return _error(f"Maximum {max_batch} variants per request; use /api/v1/annotate for VCFs", 400)

    queries = [e.variant for e in entries if e.variant]
    rows = _db().lookup_chunked(queries) if queries else []
    results = [_batch_row(r) for r in rows]
    skipped = [{"line": e.line_no, "input": f"{e.chrom} {e.pos} {e.ref} {e.alt}", "status": e.status,
                "message": e.message} for e in entries if not e.variant]
    found = sum(1 for r in results if r["found"])
    return jsonify({
        "status": "success", "genome_build": GENOME_BUILD, "count": len(results),
        "summary": {
            "submitted": len(entries), "queried": len(results), "found": found,
            "not_found": len(results) - found, "skipped": len(skipped),
            "damaging": sum(1 for r in results if r[PRED_FIELD] == "D"),
            "tolerated": sum(1 for r in results if r[PRED_FIELD] == "T"),
        },
        "not_found_variants": [{k: r[k] for k in ("chr", "pos", "ref", "alt")} for r in results if not r["found"]],
        "skipped": skipped,
        "results": results,
    })


@api_bp.post("/annotate")
@limiter.limit(_limit)
def api_annotate():
    upload = request.files.get("file") or request.files.get("vcf") or request.files.get("variants_file")
    if upload and upload.filename:
        text = decode_upload(upload.read())
        name = Path(upload.filename).name
    elif request.form.get("text", "").strip():
        text, name = request.form["text"], "variants.txt"
    elif request.data:
        text, name = decode_upload(request.data), "variants.vcf"
    else:
        return _error("Upload a VCF (form field 'file') or send text in 'text'", 400)

    fmt = (request.form.get("format") or request.args.get("format") or "json").lower()
    if fmt not in {"json", "vcf", "tsv"}:
        return _error("format must be json, vcf or tsv", 400)
    try:
        started = time.time()
        result = annotate(text, _db(), max_alleles=current_app.config["MAX_VCF_ALLELES"])
    except ValueError as exc:
        return _error(str(exc), 400)

    stem = name.split(".")[0] or "variants"
    if fmt == "vcf":
        return Response(to_vcf(result), mimetype="text/plain; charset=utf-8",
                        headers={"Content-Disposition": f"attachment; filename={stem}.annotatemissense.vcf"})
    if fmt == "tsv":
        return Response(to_tsv(result), mimetype="text/tab-separated-values; charset=utf-8",
                        headers={"Content-Disposition": f"attachment; filename={stem}.annotatemissense.tsv"})
    payload = to_json(result, limit=int(request.args.get("limit", "5000")))
    payload["seconds"] = round(time.time() - started, 3)
    return jsonify(payload)


@api_bp.get("/gene/<gene_name>")
@limiter.limit(_limit)
def api_gene(gene_name: str):
    try:
        limit = max(1, min(500, int(request.args.get("limit", "100"))))
    except ValueError:
        return _error("limit must be an integer", 400)
    db = _db()
    if not db.has_annotations:
        return _error("Gene search needs the annotation database (variants.duckdb)", 501)
    results = [_batch_row(r) for r in db.gene_variants(gene_name, limit)]
    return jsonify({"status": "success", "gene": gene_name, "count": len(results), "results": results})


@api_bp.post("/compare")
@limiter.limit(_limit)
def compare_variants():
    payload = request.get_json(silent=True) or {}
    raw = payload.get("variants", request.form.get("variants"))
    if isinstance(raw, list):
        text = "\n".join(" ".join(str(v.get(k, "")) for k in ("chr", "pos", "ref", "alt")) if isinstance(v, dict)
                         else str(v) for v in raw)
    elif isinstance(raw, str):
        try:
            decoded = json.loads(raw)
            text = "\n".join(map(str, decoded)) if isinstance(decoded, list) else raw
        except json.JSONDecodeError:
            text = raw
    else:
        return _error("Provide 'variants' (2-4 variants)", 400)
    try:
        variants = [e.variant for e in parse(text).entries if e.variant]
    except ValueError as exc:
        return _error(str(exc), 400)
    if not 2 <= len(variants) <= 4:
        return _error("Provide between 2 and 4 single-nucleotide variants", 400)

    rows = _db().lookup(variants, ("Gene.refGene", "CLNSIG", "CLNDN", "CLNREVSTAT", "am_pathogenicity",
                                   "CADD_phred", "REVEL_score", "MetaRNN_score", "BayesDel_addAF_score",
                                   "SIFT_score", "Polyphen2_HDIV_score"))
    out, best_key, best = [], None, float("-inf")
    for r in rows:
        key = f"{r['query_chr']}:{r['query_pos']} {r['query_ref']}>{r['query_alt']}"
        score = r.get(SCORE_FIELD)
        if score is not None and score > best:
            best, best_key = score, key
        votes = [(r.get(PRED_FIELD) == "D") if r.get(PRED_FIELD) else None]
        for value, threshold in ((r.get("am_pathogenicity"), 0.5), (r.get("REVEL_score"), 0.5),
                                 (r.get("MetaRNN_score"), 0.5), (r.get("CADD_phred"), 20.0)):
            votes.append(None if value is None else float(value) >= threshold)
        cast = [v for v in votes if v is not None]
        out.append({
            "chr": r["query_chr"], "pos": r["query_pos"], "ref": r["query_ref"], "alt": r["query_alt"],
            "found": r["found"], "gene": r.get("Gene.refGene"),
            "scores": {SCORE_FIELD: score, PRED_FIELD: r.get(PRED_FIELD), "AlphaMissense": r.get("am_pathogenicity"),
                       "CADD_phred": r.get("CADD_phred"), "REVEL_score": r.get("REVEL_score"),
                       "MetaRNN_score": r.get("MetaRNN_score"), "BayesDel_addAF_score": r.get("BayesDel_addAF_score"),
                       "SIFT_score": r.get("SIFT_score"), "Polyphen2_HDIV_score": r.get("Polyphen2_HDIV_score")},
            "clinical": {"CLNSIG": r.get("CLNSIG"), "CLNDN": r.get("CLNDN"), "CLNREVSTAT": r.get("CLNREVSTAT")},
            "tool_agreement": {"damaging_votes": sum(cast), "tools_with_scores": len(cast)},
        })
    return jsonify({"status": "success", "variants": out,
                    "comparison": {"highest_annotatemissense_score": best_key,
                                   "tools": "AnnotateMissense, AlphaMissense>=0.5, REVEL>=0.5, MetaRNN>=0.5, CADD>=20"}})


@api_bp.get("/stats")
@limiter.limit(_limit)
def get_stats():
    cached = _CACHE.get("stats")
    if cached and time.time() < cached[0]:
        return jsonify(cached[1])
    payload = {"status": "success", "genome_build": GENOME_BUILD, **_db().stats()}
    _CACHE["stats"] = (time.time() + _CACHE_TTL, payload)
    return jsonify(payload)


@api_bp.get("/schema")
@limiter.limit(_limit)
def get_schema():
    return jsonify({"status": "success", "tables": _db().schema(),
                    "legacy_names": {"UQ_score": SCORE_FIELD, "UQ_sd": CONFIDENCE_FIELD, "UQ_pred": PRED_FIELD}})


@api_bp.get("/sample/<fmt>")
def get_sample(fmt: str):
    files = {"txt": "sample_single.txt", "tsv": "sample_batch.tsv", "vcf": "sample_vcf.vcf"}
    filename = files.get(fmt.lower())
    if not filename:
        return _error("Unsupported sample format. Use txt, tsv, or vcf.", 400)
    samples = Path(current_app.root_path) / "static" / "samples"
    return send_from_directory(samples, filename, as_attachment=True, download_name=filename)
