import os
from pathlib import Path

import duckdb
import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def annotation_db(tmp_path_factory) -> str:
    """A 146-row extract of the released 'variants' table (chr22), as a DuckDB file."""
    path = tmp_path_factory.mktemp("db") / "variants_mini.duckdb"
    con = duckdb.connect(str(path))
    con.execute(f"CREATE TABLE variants AS SELECT * FROM read_parquet('{(DATA / 'mini_variants.parquet').as_posix()}')")
    con.close()
    return str(path)


@pytest.fixture(scope="session")
def score_index(tmp_path_factory) -> str:
    from annotatemissense.db import build_score_index

    path = tmp_path_factory.mktemp("scores") / "annotatemissense_scores.duckdb"
    build_score_index(str(DATA / "mini_scores.tsv.gz"), str(path))
    return str(path)


@pytest.fixture()
def client(annotation_db, monkeypatch):
    monkeypatch.setenv("ANNOTATEMISSENSE_DB", annotation_db)
    monkeypatch.delenv("ANNOTATEMISSENSE_SCORES", raising=False)
    monkeypatch.setenv("RATE_LIMIT", "10000/hour")
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    return app.test_client()
