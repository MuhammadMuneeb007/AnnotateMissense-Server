"""AnnotateMissense web server and REST API (Flask)."""

import os

from flask import Flask, jsonify, request
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from annotatemissense import GENOME_BUILD, ZENODO_DOI, __version__
from annotatemissense.db import AnnotateMissenseDB

cors = CORS(resources={r"/api/*": {"origins": "*"}})
limiter = Limiter(key_func=get_remote_address, storage_uri="memory://")

REPO_URL = "https://github.com/MuhammadMuneeb007/AnnotateMissense"
SERVER_REPO_URL = "https://github.com/MuhammadMuneeb007/AnnotateMissense-Server"
DOCS_URL = "https://muhammadmuneeb007.github.io/AnnotateMissense/"


def _default_path(root: str, *candidates: str):
    for name in candidates:
        path = os.path.join(root, "data", name)
        if os.path.isfile(path):
            return path
    return None


def _coverage(db: AnnotateMissenseDB) -> dict:
    """Summarise what this instance can answer (shown in the UI and /health)."""
    stats = db.stats()
    ann = stats.get("annotation_database") or {}
    idx = stats.get("score_index") or {}
    chroms = sorted(set(ann.get("chromosomes", {})) | set(idx.get("chromosomes", {})),
                    key=lambda c: (len(c), c))
    genome_wide = len(chroms) >= 24
    return {
        "chromosomes": chroms,
        "genome_wide": genome_wide,
        "annotation_rows": ann.get("rows"),
        "annotation_columns": ann.get("columns"),
        "score_index_alleles": idx.get("alleles"),
        "label": ("all chromosomes" if genome_wide else ", ".join(c.replace("chr", "chr ") for c in chroms)),
    }


def create_app() -> Flask:
    app = Flask(__name__)
    root = os.path.abspath(os.path.join(app.root_path, ".."))

    annotations = os.getenv("ANNOTATEMISSENSE_DB") or _default_path(root, "variants.duckdb", "variants_dev.duckdb")
    scores = os.getenv("ANNOTATEMISSENSE_SCORES") or _default_path(root, "annotatemissense_scores.duckdb")
    db = AnnotateMissenseDB(annotations, scores)

    app.config.update(
        VARIANT_DB=db,
        VERSION=__version__,
        MAX_CONTENT_LENGTH=int(os.getenv("MAX_UPLOAD_MB", "50")) * 1024 * 1024,
        MAX_BATCH=int(os.getenv("MAX_BATCH_VARIANTS", "1000")),
        MAX_VCF_ALLELES=int(os.getenv("MAX_VCF_ALLELES", "10000")),
        RATE_LIMIT=os.getenv("RATE_LIMIT", "300/hour"),
        COVERAGE=_coverage(db),
        INSTANCE_NOTE=os.getenv("INSTANCE_NOTE", ""),
    )

    if os.getenv("TRUST_PROXY", "0") == "1":
        from werkzeug.middleware.proxy_fix import ProxyFix

        app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1)

    cors.init_app(app)
    limiter.init_app(app)

    from .api import api_bp
    from .routes import web_bp

    app.register_blueprint(web_bp)
    app.register_blueprint(api_bp)

    @app.context_processor
    def inject_globals():
        return {
            "coverage": app.config["COVERAGE"],
            "instance_note": app.config["INSTANCE_NOTE"],
            "app_version": __version__,
            "genome_build": GENOME_BUILD,
            "zenodo_doi": ZENODO_DOI,
            "repo_url": REPO_URL,
            "server_repo_url": SERVER_REPO_URL,
            "docs_url": DOCS_URL,
            "max_batch": app.config["MAX_BATCH"],
            "max_vcf_alleles": app.config["MAX_VCF_ALLELES"],
        }

    @app.after_request
    def apply_security_headers(response):
        response.headers.pop("Set-Cookie", None)
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        return response

    def _api_error(message: str, code: int):
        return jsonify({"status": "error", "message": message}), code

    @app.errorhandler(404)
    def handle_404(error):
        if request.path.startswith("/api/"):
            return _api_error("Resource not found", 404)
        return error

    @app.errorhandler(405)
    def handle_405(error):
        if request.path.startswith("/api/"):
            return _api_error("Method not allowed", 405)
        return error

    @app.errorhandler(413)
    def handle_413(error):
        return _api_error(f"Upload too large (limit {os.getenv('MAX_UPLOAD_MB', '50')} MB)", 413)

    @app.errorhandler(429)
    def handle_429(error):
        return _api_error(f"Rate limit exceeded: {error.description}", 429)

    @app.errorhandler(Exception)
    def handle_exception(error):
        from werkzeug.exceptions import HTTPException

        if isinstance(error, HTTPException):
            if request.path.startswith("/api/"):
                return _api_error(error.description or error.name, error.code or 500)
            return error
        if request.path.startswith("/api/"):
            app.logger.exception("Unhandled API error")
            return _api_error("Internal server error", 500)
        raise error

    return app
