#!/usr/bin/env python3
"""Build the small demo dataset used by the hosted (Vercel) instance.

From a chromosome-22 extract of the released database (variants_dev.duckdb) it writes:

  demo/annotatemissense_demo_scores.duckdb       score index for every chr22 allele (table `scores`)
  demo/annotatemissense_demo_annotations.duckdb  all 337 columns for a panel of chr22 disease genes
                                                 (table `variants`, same layout as the release)

    python scripts/build_demo_data.py data/variants_dev.duckdb --out demo
"""

import argparse
import os
from pathlib import Path

import duckdb

PANEL = ["MYH9", "NF2", "CHEK2", "SMARCB1", "COMT", "LZTR1", "ADSL", "UBE2L3", "TBX1", "SOX10",
         "ARSA", "TCN2", "PDGFB", "CRYBB2", "SLC5A1"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", help="variants_dev.duckdb (chr22 extract of variants.duckdb)")
    parser.add_argument("--out", default="demo")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scores_path = out / "annotatemissense_demo_scores.duckdb"
    ann_path = out / "annotatemissense_demo_annotations.duckdb"
    for p in (scores_path, ann_path):
        if p.exists():
            p.unlink()

    src = Path(args.source).resolve().as_posix()

    con = duckdb.connect(str(scores_path))
    con.execute(f"ATTACH '{src}' AS src (READ_ONLY)")
    con.execute("""
        CREATE TABLE scores AS
        SELECT Chr AS chrom, CAST(Start AS INTEGER) AS pos, Ref AS ref, Alt AS alt,
               CAST(any_value(UQ_score) AS FLOAT) AS score, CAST(any_value(UQ_sd) AS FLOAT) AS confidence,
               any_value(UQ_pred) AS pred, CAST(COUNT(*) AS USMALLINT) AS source_rows
        FROM src.variants GROUP BY ALL ORDER BY chrom, pos, ref, alt
    """)
    n = con.execute("SELECT COUNT(*), COUNT(score) FROM scores").fetchone()
    con.execute("CREATE TABLE meta (key VARCHAR, value VARCHAR)")
    con.executemany("INSERT INTO meta VALUES (?, ?)", [
        ("source_file", os.path.basename(args.source)), ("source_columns", "UQ_score, UQ_sd, UQ_pred"),
        ("alleles", str(n[0])), ("scored_alleles", str(n[1])), ("genome_build", "GRCh38"),
    ])
    con.close()

    con = duckdb.connect(str(ann_path))
    con.execute(f"ATTACH '{src}' AS src (READ_ONLY)")
    con.execute("""
        CREATE TABLE variants AS
        SELECT * FROM src.variants
        WHERE list_has_any(string_split(regexp_replace("Gene.refGene", ';', ',', 'g'), ','), ?)
        ORDER BY Chr, Start, rowid
    """, [PANEL])
    rows = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    con.close()

    for p in (scores_path, ann_path):
        print(f"{p}: {p.stat().st_size / 1e6:.1f} MB")
    print(f"scores: {n[0]:,} alleles ({n[1]:,} scored); annotations: {rows:,} rows for {len(PANEL)} genes")


if __name__ == "__main__":
    main()
