import pytest

from clade.provenance.manifest import ProvenanceError, register_file, validate_manifest


def test_clean_pass(tmp_path):
    data_file = tmp_path / "data.tsv"
    data_file.write_text("sample\tvalue\nA\t1\n")
    manifest_path = tmp_path / "manifest.yaml"

    register_file(manifest_path, data_file, generating_script="derive_phenotype.py")

    # Should not raise.
    validate_manifest(manifest_path, [data_file])


def test_missing_manifest_entry_raises(tmp_path):
    data_file = tmp_path / "unregistered.tsv"
    data_file.write_text("sample\tvalue\nA\t1\n")
    manifest_path = tmp_path / "manifest.yaml"
    manifest_path.write_text("files: {}\n")

    with pytest.raises(ProvenanceError, match="no manifest entry"):
        validate_manifest(manifest_path, [data_file])


def test_checksum_mismatch_raises(tmp_path):
    data_file = tmp_path / "phenotype.tsv"
    data_file.write_text("sample\tphenotype\nA\t1\n")
    manifest_path = tmp_path / "manifest.yaml"

    register_file(manifest_path, data_file, generating_script="derive_phenotype.py")

    # Mutate the file after registration -- the exact failure mode a
    # mislabeled/swapped phenotype file would trigger.
    data_file.write_text("sample\tphenotype\nA\t0\n")

    with pytest.raises(ProvenanceError, match="checksum mismatch"):
        validate_manifest(manifest_path, [data_file])


def test_entry_without_provenance_declaration_raises(tmp_path):
    data_file = tmp_path / "no_provenance.tsv"
    data_file.write_text("sample\tvalue\nA\t1\n")
    manifest_path = tmp_path / "manifest.yaml"

    with pytest.raises(ProvenanceError, match="must supply generating_script or external_source"):
        register_file(manifest_path, data_file)


def test_register_file_requires_existing_file(tmp_path):
    manifest_path = tmp_path / "manifest.yaml"
    with pytest.raises(FileNotFoundError):
        register_file(manifest_path, tmp_path / "does_not_exist.tsv", generating_script="x.py")
