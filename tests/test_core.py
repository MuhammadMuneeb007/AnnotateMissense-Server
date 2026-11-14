import pytest

from annotatemissense.db import (AnnotateMissenseDB, normalize_chr, normalize_variant, parse_variant_string,
                                 trim_alleles)
from annotatemissense.vcf import annotate, parse, to_tsv, to_vcf

MYH9_R702C = ("chr22", 36305985, "G", "A")
MULTI_ROW = ("chr22", 21549460, "C", "A")


def test_normalisation():
    assert normalize_chr("22") == normalize_chr("chr22") == normalize_chr("CHR22") == "chr22"
    assert normalize_chr("MT") == normalize_chr("chrM") == "chrM"
    assert normalize_variant("22", "36305985", "g", "a") == MYH9_R702C
    assert parse_variant_string("22:36305985:G:A") == parse_variant_string("chr22-36305985-G-A") == MYH9_R702C
    assert trim_alleles(36305984, "CG", "CA") == (36305985, "G", "A")
    with pytest.raises(ValueError):
        normalize_chr("23")
    with pytest.raises(ValueError):
        normalize_variant("1", "x", "A", "G")


def test_lookup_one_record_per_allele(annotation_db):
    db = AnnotateMissenseDB(annotation_db)
    rows = db.lookup([MYH9_R702C, MULTI_ROW, ("chr22", 1, "A", "G")])
    assert [r["found"] for r in rows] == [True, True, False]
    assert rows[0]["annotatemissense_pred"] == "D" and rows[0]["annotatemissense_score"] > 0.99
    assert rows[1]["annotation_rows"] == 2  # two ESM protein joins, one allele
    assert rows[2]["annotatemissense_score"] is None


def test_transcript_context(annotation_db):
    ctx = AnnotateMissenseDB(annotation_db).transcript_context(MULTI_ROW)
    assert ctx["annotation_rows"] == 2
    assert ctx["annotatemissense_score_is_per_allele"] is True
    assert ctx["dbnsfp"]["genes"] == ["UBE2L3"]
    assert ctx["alphamissense"]["transcript"].startswith("ENST")
    assert len(ctx["esm1v"]) == 2


def test_audit_and_schema(annotation_db):
    db = AnnotateMissenseDB(annotation_db)
    audit = db.audit_alleles()
    assert audit["rows"] == 146 and audit["alleles"] == 134
    assert audit["alleles_with_discordant_scores"] == 0
    schema = db.schema()["variants"]
    assert len(schema) == 337
    assert {c["name"] for c in schema} >= {"Chr", "Start", "Ref", "Alt", "UQ_score", "am_pathogenicity"}


def test_score_index(score_index):
    db = AnnotateMissenseDB(scores_path=score_index)
    stats = db.stats()["score_index"]
    assert stats["alleles"] < 70  # duplicates collapsed
    assert stats["scored_alleles"] == stats["alleles"] - 4  # 'U' rows -> NULL score, not 0
    row = db.lookup([("chr1", 65565, "A", "T")])[0]
    assert row["found"] and row["annotatemissense_pred"] == "D" and abs(row["annotatemissense_score"] - 0.814325) < 1e-5
    unscored = db.lookup([("chr1", 941172, "G", "A")])[0]
    assert unscored["found"] and unscored["annotatemissense_score"] is None


def test_scores_and_annotations_together(annotation_db, score_index):
    db = AnnotateMissenseDB(annotation_db, score_index)
    rows = db.lookup([MYH9_R702C, ("chr1", 65565, "A", "T")])
    assert rows[0]["found"] and rows[0]["Gene.refGene"] == "MYH9"
    assert rows[1]["found"] and rows[1]["Gene.refGene"] is None and rows[1]["in_score_index"]


VCF = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr22\t21549460\t.\tC\tA,G\t.\tPASS\tDP=10
chr22\t36305984\t.\tCG\tCA\t.\tPASS\t.
chr22\t36305985\t.\tGG\tG\t.\tPASS\t.
chr22\t10000000\t.\tA\tG\t.\tPASS\t.
chr99\t1\t.\tA\tG\t.\tPASS\t.
"""


def test_vcf_annotation(annotation_db):
    result = annotate(VCF, AnnotateMissenseDB(annotation_db))
    assert [e.status for e in result.entries] == ["found", "found", "found", "not_snv", "not_in_release", "invalid"]
    summary = result.summary()
    assert summary["records"] == 5 and summary["alt_alleles"] == 6
    out = to_vcf(result)
    assert "##INFO=<ID=AnnotateMissense_score,Number=A" in out
    multi = [l for l in out.splitlines() if l.startswith("chr22\t21549460")][0]
    info = dict(kv.split("=") for kv in multi.split("\t")[7].split(";"))
    assert info["DP"] == "10" and len(info["AnnotateMissense_score"].split(",")) == 2
    assert info["AnnotateMissense_status"] == "found,found"
    assert to_tsv(result).count("\n") == 7


def test_grch37_warning():
    result = parse("##fileformat=VCFv4.1\n##reference=GRCh37/hg19\n#CHROM\tPOS\tID\tREF\tALT\n22\t1\t.\tA\tG\n")
    assert result.warnings and "GRCh37" in result.warnings[0]


def test_variant_list_and_limit():
    result = parse("chrom pos ref alt\n22 36305985 G A\n22,19963748,G,A\n")
    assert [e.variant for e in result.entries] == [MYH9_R702C, ("chr22", 19963748, "G", "A")]
    with pytest.raises(ValueError):
        parse("22 1 A G\n22 2 A G\n", max_alleles=1)


def test_large_lookup_path_matches_inline(annotation_db):
    db = AnnotateMissenseDB(annotation_db)
    variants = [MYH9_R702C, MULTI_ROW, ("chr22", 1, "A", "G"), ("chr22", 19963748, "G", "A")]
    inline = db.lookup(variants)
    via_temp_table = db.lookup_chunked(variants, chunk_size=2)
    assert [r["idx"] for r in via_temp_table] == [0, 1, 2, 3]
    assert [(r["found"], r["annotatemissense_score"]) for r in inline] == \
           [(r["found"], r["annotatemissense_score"]) for r in via_temp_table]
