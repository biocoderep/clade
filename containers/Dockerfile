# Container for the CLADE package and the pipeline's Python steps.
#
# Pinned to Python 3.10 for one specific reason: the optional `firth` extra
# depends on firthlogist, which uses a private scikit-learn API (_validate_data)
# removed in scikit-learn 1.6. That constrains scikit-learn to <1.6, which in
# turn constrains Python. Pinning here rather than documenting it means the
# constraint cannot be violated by accident.
FROM python:3.10-slim-bookworm

LABEL org.opencontainers.image.title="clade" \
      org.opencontainers.image.description="Lineage-aware validation of candidate compensatory mutations" \
      org.opencontainers.image.source="https://github.com/biocoderep/clade" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONNOUSERSITE=1

# procps supplies `ps`, which Nextflow calls to collect per-task resource
# metrics. Without it every task logs a warning and the trace file has empty
# memory columns.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps bash ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/clade
COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir ".[firth]" pyyaml \
    && pip cache purge

# The pipeline's helper scripts, on PATH so modules call them by name.
COPY workflows/nextflow/bin/ /usr/local/bin/
RUN chmod +x /usr/local/bin/*.py

RUN clade --version && python -c "import clade, dendropy, statsmodels, pandas, numpy, scipy"

CMD ["clade", "--help"]
