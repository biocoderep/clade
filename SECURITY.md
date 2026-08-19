# Security Policy

CLADE is a research-analysis tool, not a network service — the realistic security surface is small (it reads local files and calls a handful of well-established scientific Python libraries). That said:

## Reporting a vulnerability

If you find a security issue (e.g. unsafe deserialization, a dependency with a known CVE that CLADE pulls in), please email `bioinfocode4@gmail.com` rather than opening a public issue, so there's time to address it before disclosure.

## Scope notes

- CLADE does not execute arbitrary user-supplied code paths beyond loading the CSV/TSV/Newick files a user explicitly points it at.
- CLADE never phones home, uploads data, or makes network calls itself — genome retrieval (`prefetch`/`fasterq-dump`) is an external, separate step in the genome-processing pipeline, not part of the `clade` package.
