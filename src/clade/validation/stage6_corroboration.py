"""
Stage 6 — External corroboration.

Unlike Stages 1-5, this stage is not automated in CLADE, and we do not
pretend otherwise: it is a literature search and, where the original locus
annotation is generic, a homology search (e.g. BLASTp) to resolve the
candidate's specific identity first. In the CLADE case study, this was run
manually (BLASTp via NCBI's web interface, PubMed search by hand) — no
script produced these results, and none should claim to.

This module exists only to give that manual process a consistent, explicit
place to be recorded, so Stage 6 findings are structured data (queryable,
diffable) rather than prose scattered across a manuscript.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CorroborationRecord:
    """A structured record of one candidate's Stage 6 finding.

    `mechanistic_support`: True/False/None (None = search performed, no
    clear answer either way — distinct from "search not yet performed",
    which should simply not have a CorroborationRecord at all).
    """

    candidate: str
    resolved_identity: str | None = None
    identity_resolution_method: str | None = None  # e.g. "BLASTp vs NCBI nr, 2026-08-12"
    literature_support: str | None = None  # e.g. a PMID, or "none found"
    mechanistic_support: bool | None = None
    notes: str = ""
    sources: list[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.mechanistic_support is True:
            return f"{self.candidate}: mechanistic support found ({self.literature_support})"
        if self.mechanistic_support is False:
            return f"{self.candidate}: no mechanistic support found despite search ({self.identity_resolution_method or 'no identity resolution attempted'})"
        return f"{self.candidate}: Stage 6 inconclusive"
