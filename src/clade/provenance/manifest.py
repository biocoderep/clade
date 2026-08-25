"""
Provenance manifest: structural enforcement for the practice already stated
in docs/reproducibility/phenotype_correction.md and
docs/reproducibility/amr_panel_provenance.md ("derive files programmatically;
treat any file without a documented generating script as unconfirmed"), but
not previously enforced by any tooling — both of those incidents were caught
by hand. `validate_manifest` raises rather than warns, by design: a file
used under a wrong label fails the checksum check; a file with no recorded
generating script or external source fails the provenance check. Neither
class of error can pass silently through a call to `validate_manifest`.

This module registers files going forward. It does not retroactively cover
every file already in this repository — see the manifest file itself
(`provenance_manifest.yaml`, repo root) for exactly which files have been
registered so far, and treat every unregistered file exactly as
`docs/reproducibility/*.md` already recommend: unconfirmed until checked.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import yaml


class ProvenanceError(Exception):
    """Raised when a file fails checksum verification or has no (or an
    incomplete) manifest entry. Callers should let this propagate, not
    catch-and-continue -- a passing analysis step is exactly the thing
    this module exists to gate."""


@dataclass
class ManifestEntry:
    path: str
    sha256: str
    generating_script: str | None = None
    external_source: str | None = None
    note: str | None = None

    def has_provenance(self) -> bool:
        return bool(self.generating_script or self.external_source)


def _sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest(manifest_path: str | Path) -> dict[str, ManifestEntry]:
    """Load a manifest YAML file into {path: ManifestEntry}. A missing
    manifest file loads as empty, not an error -- validate_manifest is
    where "no entry for this file" becomes a real failure."""
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        return {}
    with open(manifest_path) as f:
        raw = yaml.safe_load(f) or {}
    entries = {}
    for path, fields in (raw.get("files") or {}).items():
        entries[path] = ManifestEntry(
            path=path,
            sha256=fields.get("sha256", ""),
            generating_script=fields.get("generating_script"),
            external_source=fields.get("external_source"),
            note=fields.get("note"),
        )
    return entries


def register_file(
    manifest_path: str | Path,
    file_path: str | Path,
    generating_script: str | None = None,
    external_source: str | None = None,
    note: str | None = None,
) -> ManifestEntry:
    """Compute `file_path`'s real sha256 now and record it in the manifest,
    along with a mandatory provenance declaration.

    Requires `generating_script` or `external_source` (at least one) --
    register_file refuses a file with neither, for the same reason
    validate_manifest would later reject it: an undocumented file is
    exactly the failure mode this module exists to catch, not defer.
    """
    if generating_script is None and external_source is None:
        raise ProvenanceError(
            f"register_file({file_path}): must supply generating_script or "
            f"external_source. A file registered with neither would still "
            f"fail validate_manifest's provenance check later -- refusing "
            f"to create that state now."
        )

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    manifest_path = Path(manifest_path)
    raw: dict = {}
    if manifest_path.exists():
        with open(manifest_path) as f:
            raw = yaml.safe_load(f) or {}
    raw.setdefault("files", {})

    key = str(file_path)
    entry = ManifestEntry(
        path=key,
        sha256=_sha256_of(file_path),
        generating_script=generating_script,
        external_source=external_source,
        note=note,
    )
    raw["files"][key] = {
        "sha256": entry.sha256,
        "generating_script": entry.generating_script,
        "external_source": entry.external_source,
        "note": entry.note,
    }

    with open(manifest_path, "w") as f:
        yaml.safe_dump(raw, f, sort_keys=True)

    return entry


def validate_manifest(manifest_path: str | Path, files: list[str | Path]) -> None:
    """Verify every path in `files` against the manifest at `manifest_path`.

    Raises `ProvenanceError`, naming the exact file and exact reason, on
    the first problem found, in this order:
      1. no manifest entry at all for the file
      2. an entry exists but declares neither generating_script nor
         external_source
      3. the file's current sha256 does not match the manifest's recorded
         sha256 (the file changed, or a different file was substituted,
         since it was registered)

    Returns None (no error) if every file passes all three checks. Intended
    to be called at the top of any analysis step, before that step reads
    its inputs -- see module docstring for why this needs to raise, not warn.
    """
    entries = load_manifest(manifest_path)

    for f in files:
        f = str(f)
        if f not in entries:
            raise ProvenanceError(
                f"{f}: no manifest entry. Register it with register_file() "
                f"before using it in any analysis step."
            )
        entry = entries[f]
        if not entry.has_provenance():
            raise ProvenanceError(
                f"{f}: manifest entry exists but declares no generating_script "
                f"or external_source -- undocumented provenance, treat as "
                f"unconfirmed (see docs/reproducibility/phenotype_correction.md)."
            )
        actual = _sha256_of(Path(f))
        if actual != entry.sha256:
            raise ProvenanceError(
                f"{f}: checksum mismatch. Manifest recorded {entry.sha256}, "
                f"file on disk is now {actual}. This file changed, or a "
                f"different file was substituted, since it was registered -- "
                f"this is exactly the failure mode that let a blaOXA-66 file "
                f"get used under the blaOXA-23 label "
                f"(docs/reproducibility/phenotype_correction.md). Do not "
                f"proceed without deliberately re-registering."
            )
