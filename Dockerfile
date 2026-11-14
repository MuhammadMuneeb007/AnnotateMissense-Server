FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /srv/annotatemissense

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY annotatemissense/ ./annotatemissense/
COPY app/ ./app/
COPY wsgi.py pyproject.toml README.md ./

# Mount the data directory at /data (variants.duckdb and/or annotatemissense_scores.duckdb)
ENV ANNOTATEMISSENSE_DB=/data/variants.duckdb
EXPOSE 8000
CMD ["gunicorn", "--bind", "0.0.0.0:8000", "--workers", "2", "--threads", "8", "--timeout", "300", "--access-logfile", "-", "wsgi:app"]
