// Collect every process's versions.yml into one file.
//
// This is the artefact a reader needs in order to reproduce a result: not
// "we used FastQC" but the exact version of every tool that touched the data
// in THIS run. Without it a pipeline is reproducible only until someone
// changes a container tag.
//
// Deliberately parses the simple two-level YAML by hand rather than importing
// PyYAML: the stock python container has no third-party packages, and adding a
// dependency here would mean this step could fail for reasons unrelated to the
// science it is recording.
process COLLATE_VERSIONS {
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "quay.io/biocontainers/python:3.10"

    input:
    path versions

    output:
    path "software_versions.yml", emit: yml

    script:
    """
    #!/usr/bin/env python3
    import glob, platform

    merged = {}
    for fn in sorted(glob.glob("*.yml")) + sorted(glob.glob("*.yaml")):
        section = None
        for raw in open(fn):
            line = raw.rstrip("\\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if not line.startswith((" ", "\\t")):
                # strip the trailing colon, NOT split on ":" -- process
                # names are themselves colon-separated (PER_SAMPLE:FASTP), and
                # splitting collapses every process in a subworkflow into one.
                section = line.strip().rstrip(":").strip().strip('"')
                merged.setdefault(section, {})
            elif section and ":" in line:
                k, _, v = line.strip().partition(":")
                merged[section][k.strip().strip('"')] = v.strip().strip('"')

    merged["Workflow"] = {
        "biocoderep/clade": "${workflow.manifest.version}",
        "Nextflow": "${workflow.nextflow.version}",
        "Python": platform.python_version(),
    }

    with open("software_versions.yml", "w") as out:
        for section in sorted(merged):
            out.write('"%s":\\n' % section)
            for tool in sorted(merged[section]):
                out.write("    %s: %s\\n" % (tool, merged[section][tool]))
    """

    stub:
    """
    printf '"Workflow":\\n    biocoderep/clade: ${workflow.manifest.version}\\n' > software_versions.yml
    """
}
