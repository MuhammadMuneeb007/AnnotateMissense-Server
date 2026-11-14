"""Command-line interface: ``annotatemissense <command> ...``

Commands
--------
lookup              Look up one or more alleles (CHR:POS:REF:ALT).
annotate-vcf        Annotate a VCF / variant list offline (no size limit, nothing uploaded).
schema              Print the table schema of the database(s).
stats               Print row / allele / prediction counts.
audit               Allele-multiplicity audit (rows per allele, discordant scores).
build-score-index   Build the compact score index from UQ_BioSig_model_Final.tsv.gz.
serve               Run the web server and REST API.

Databases are given with ``--db`` (annotation database, variants.duckdb) and/or
``--scores`` (score index), or the ANNOTATEMISSENSE_DB / ANNOTATEMISSENSE_SCORES
environment variables.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from typing import List, Optional

from . import CONFIDENCE_FIELD, PRED_FIELD, SCORE_FIELD, __version__
from .db import AnnotateMissenseDB, build_score_index, parse_variant_string


def _open_db(args: argparse.Namespace) -> AnnotateMissenseDB:
    db_path = args.db or os.environ.get("ANNOTATEMISSENSE_DB")
    scores = args.scores or os.environ.get("ANNOTATEMISSENSE_SCORES")
    if not db_path and not scores:
        sys.exit("error: give --db variants.duckdb and/or --scores annotatemissense_scores.duckdb")
    return AnnotateMissenseDB(db_path, scores)


def _read_input(path: str) -> str:
    from .vcf import decode_upload

    if path == "-":
        return decode_upload(sys.stdin.buffer.read())
    with open(path, "rb") as handle:
        return decode_upload(handle.read())


def cmd_lookup(args: argparse.Namespace) -> None:
    db = _open_db(args)
    raw: List[str] = list(args.variants)
    if args.file:
        raw += [line.strip() for line in open(args.file, encoding="utf-8") if line.strip() and not line.startswith("#")]
    if not raw:
        sys.exit("error: give variants as CHR:POS:REF:ALT or --file")
    variants = [parse_variant_string(v) for v in raw]
    rows = db.lookup_chunked(variants)
    if args.json:
        print(json.dumps(rows, indent=2, default=str))
        return
    cols = ["query_chr", "query_pos", "query_ref", "query_alt", "found", SCORE_FIELD, CONFIDENCE_FIELD,
            PRED_FIELD, "Gene.refGene", "AAChange.refGene", "CLNSIG", "am_pathogenicity", "annotation_rows"]
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(cols)
    for row in rows:
        writer.writerow(["" if row.get(c) is None else row.get(c) for c in cols])


def cmd_annotate_vcf(args: argparse.Namespace) -> None:
    from .vcf import annotate, to_json, to_tsv, to_vcf

    db = _open_db(args)
    start = time.time()
    result = annotate(_read_input(args.input), db)
    fmt = args.format or ("tsv" if args.output and args.output.endswith(".tsv") else "vcf")
    text = to_vcf(result) if fmt == "vcf" else to_tsv(result) if fmt == "tsv" else json.dumps(to_json(result), indent=2, default=str)
    if args.output and args.output != "-":
        with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    summary = result.summary()
    for warning in result.warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    print(
        "AnnotateMissense: {records} records, {alt_alleles} ALT alleles -> found {found}, "
        "not_in_release {not_in_release}, not_snv {not_snv}, no_score {no_score}, invalid {invalid} "
        "(D={damaging}, T={tolerated}) in {secs:.2f}s".format(secs=time.time() - start, **summary),
        file=sys.stderr,
    )


def cmd_schema(args: argparse.Namespace) -> None:
    db = _open_db(args)
    schema = db.schema()
    if args.json:
        print(json.dumps(schema, indent=2))
        return
    for table, cols in schema.items():
        print(f"## {table} ({len(cols)} columns)")
        for c in cols:
            print(f"{c['position']:>4}  {c['name']:<48} {c['type']:<10} {c.get('source', '')}")


def cmd_stats(args: argparse.Namespace) -> None:
    print(json.dumps(_open_db(args).stats(), indent=2))


def cmd_audit(args: argparse.Namespace) -> None:
    start = time.time()
    report = _open_db(args).audit_alleles()
    report["seconds"] = round(time.time() - start, 2)
    print(json.dumps(report, indent=2))


def cmd_build(args: argparse.Namespace) -> None:
    start = time.time()
    meta = build_score_index(args.tsv, args.output, threads=args.threads)
    meta["seconds"] = f"{time.time() - start:.1f}"
    print(json.dumps(meta, indent=2))


def cmd_serve(args: argparse.Namespace) -> None:
    if args.db:
        os.environ["ANNOTATEMISSENSE_DB"] = args.db
    if args.scores:
        os.environ["ANNOTATEMISSENSE_SCORES"] = args.scores
    try:
        from app import create_app
    except ImportError as exc:  # pragma: no cover
        sys.exit(f"error: the web server needs the server extras (pip install -r requirements.txt): {exc}")
    application = create_app()
    if args.debug:
        application.run(host=args.host, port=args.port, debug=True)
        return
    from waitress import serve

    print(f"AnnotateMissense server on http://{args.host}:{args.port}", file=sys.stderr)
    serve(application, host=args.host, port=args.port, threads=args.threads, ident="AnnotateMissense")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="annotatemissense", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"annotatemissense {__version__}")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--db", help="annotation database (variants.duckdb)")
    common.add_argument("--scores", help="score index (annotatemissense_scores.duckdb)")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("lookup", parents=[common], help="look up alleles")
    p.add_argument("variants", nargs="*", help="CHR:POS:REF:ALT, e.g. 22:15528192:A:C")
    p.add_argument("--file", help="file with one variant per line")
    p.add_argument("--json", action="store_true", help="JSON output (all returned fields)")
    p.set_defaults(func=cmd_lookup)

    p = sub.add_parser("annotate-vcf", parents=[common], help="annotate a VCF or variant list offline")
    p.add_argument("input", help="VCF / VCF.gz / TSV variant list, or - for stdin")
    p.add_argument("-o", "--output", help="output path (default stdout)")
    p.add_argument("--format", choices=["vcf", "tsv", "json"], help="output format (default vcf)")
    p.set_defaults(func=cmd_annotate_vcf)

    p = sub.add_parser("schema", parents=[common], help="print table schema")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_schema)

    p = sub.add_parser("stats", parents=[common], help="print counts")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("audit", parents=[common], help="allele-multiplicity audit (full scan)")
    p.set_defaults(func=cmd_audit)

    p = sub.add_parser("build-score-index", help="build score index from the released TSV")
    p.add_argument("tsv", help="UQ_BioSig_model_Final.tsv.gz")
    p.add_argument("-o", "--output", default="annotatemissense_scores.duckdb")
    p.add_argument("--threads", type=int)
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("serve", parents=[common], help="run the web server")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--debug", action="store_true")
    p.set_defaults(func=cmd_serve)
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except (ValueError, FileNotFoundError, FileExistsError) as exc:
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    main()
