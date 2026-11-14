from flask import Blueprint, current_app, redirect, render_template, request, url_for

from annotatemissense.db import normalize_variant

from .api import _batch_row, format_single_variant

web_bp = Blueprint("web", __name__)


@web_bp.get("/")
def index():
    return render_template("index.html", page_title="Home")


@web_bp.get("/batch")
def batch():
    return redirect(url_for("web.batch_results"))


@web_bp.get("/batch_results")
def batch_results():
    return render_template("batch_results.html", page_title="Batch results")


@web_bp.get("/annotate")
def annotate_page():
    return render_template("annotate.html", page_title="Annotate VCF")


@web_bp.get("/results")
def results():
    query = {k: request.args.get(k, "") for k in ("chr", "pos", "ref", "alt")}
    payload, error_message = None, None
    if all(query.values()):
        try:
            variant = normalize_variant(query["chr"], query["pos"], query["ref"], query["alt"])
            db = current_app.config["VARIANT_DB"]
            row = db.fetch_variant(variant)
            if row:
                payload = format_single_variant(row, db.transcript_context(variant))
            else:
                error_message = (f"{variant[0]}:{variant[1]} {variant[2]}>{variant[3]} is not among the "
                                 "precomputed GRCh38 missense alleles served here. Absence is not evidence of benignity.")
        except ValueError as exc:
            error_message = str(exc)
    return render_template("results.html", page_title="Variant", query=query, payload=payload,
                           error_message=error_message)


@web_bp.get("/api_docs")
def api_docs():
    return render_template("api_docs.html", page_title="API")


@web_bp.get("/development")
def development_page():
    return render_template("development.html", page_title="About the model")


@web_bp.get("/compare")
def compare_page():
    initial = []
    for key in ("v1", "v2", "v3", "v4"):
        parts = (request.args.get(key) or "").split(":")
        if len(parts) == 4:
            initial.append(dict(zip(("chr", "pos", "ref", "alt"), parts)))
    return render_template("compare.html", page_title="Compare", initial_variants=initial)


@web_bp.get("/gene/<gene_name>")
def gene_page(gene_name: str):
    try:
        limit = max(1, min(500, int(request.args.get("limit", "200"))))
    except ValueError:
        limit = 200
    rows = [_batch_row(r) for r in current_app.config["VARIANT_DB"].gene_variants(gene_name, limit)]
    return render_template("gene.html", page_title=gene_name, gene_name=gene_name, results=rows, limit=limit)


@web_bp.get("/help")
def help_page():
    return render_template("help.html", page_title="Help")


@web_bp.get("/privacy")
def privacy():
    return render_template("privacy.html", page_title="Privacy")


@web_bp.get("/robots.txt")
def robots():
    return ("User-agent: *\nAllow: /\nDisallow: /results\nDisallow: /gene/\nDisallow: /api/\n", 200,
            {"Content-Type": "text/plain; charset=utf-8"})
