import io


def test_health(client):
    data = client.get("/api/v1/health").get_json()
    assert data["status"] == "ok" and data["genome_build"] == "GRCh38"
    assert data["coverage"]["chromosomes"] == ["chr22"]


def test_variant_get_and_post(client):
    a = client.get("/api/v1/variant/22:36305985:G:A").get_json()
    b = client.post("/api/v1/variant", json={"chrom": "chr22", "pos": 36305985, "ref": "G", "alt": "A"}).get_json()
    assert a["annotatemissense"] == b["annotatemissense"]
    assert a["annotatemissense"]["prediction"] == "D"
    assert a["gene"]["name"] == "MYH9"
    assert a["transcripts"]["annotatemissense_score_is_per_allele"] is True


def test_variant_not_found_and_bad_input(client):
    r = client.get("/api/v1/variant/22:1:A:G")
    assert r.status_code == 404 and "not evidence of benignity" in r.get_json()["detail"]
    assert client.get("/api/v1/variant/22:abc:A:G").status_code == 400


def test_variant_full(client):
    data = client.get("/api/v1/variant/full?chrom=22&pos=36305985&ref=G&alt=A").get_json()
    assert len(data["all_columns"]) == 337 and data["all_columns"]["UQ_pred"] == "D"


def test_batch(client):
    data = client.post("/api/v1/batch", data={"variants_text": "22 36305985 G A\n22 21549460 C A,G\n22 1 A G\n22 5 AT A"}).get_json()
    assert data["summary"]["submitted"] == 5
    assert data["summary"]["found"] == 3 and data["summary"]["not_found"] == 1 and data["summary"]["skipped"] == 1


def test_annotate_formats(client):
    vcf = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr22\t36305985\t.\tG\tA\t.\t.\t.\n"
    js = client.post("/api/v1/annotate?format=json", data={"file": (io.BytesIO(vcf), "x.vcf")}).get_json()
    assert js["summary"]["found"] == 1
    out = client.post("/api/v1/annotate?format=vcf", data={"file": (io.BytesIO(vcf), "x.vcf")})
    assert out.status_code == 200 and b"AnnotateMissense_pred=D" in out.data
    tsv = client.post("/api/v1/annotate?format=tsv", data={"text": "22 36305985 G A"})
    assert tsv.data.decode().splitlines()[1].split("\t")[7] == "found"


def test_gene_compare_schema_stats(client):
    gene = client.get("/api/v1/gene/MYH9?limit=10").get_json()
    assert gene["count"] == 10 and all(r["gene"] for r in gene["results"])
    cmp = client.post("/api/v1/compare", json={"variants": ["22:36305985:G:A", "22:19963748:G:A"]}).get_json()
    assert cmp["comparison"]["highest_annotatemissense_score"].startswith("chr22:36305985")
    assert len(client.get("/api/v1/schema").get_json()["tables"]["variants"]) == 337
    assert client.get("/api/v1/stats").get_json()["annotation_database"]["rows"] == 146


def test_pages_render(client):
    for path in ["/", "/annotate", "/api_docs", "/help", "/privacy", "/development", "/compare",
                 "/batch_results", "/gene/MYH9", "/results?chr=22&pos=36305985&ref=G&alt=A"]:
        r = client.get(path)
        assert r.status_code == 200, path
        assert b"AnnotateMissense" in r.data
        assert "Set-Cookie" not in r.headers
