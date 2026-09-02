#!/bin/bash
###############################################################################
# 7_Environment_Reproducibility.sh
#
# INTERACTIVE helper — run this directly on the login node (do NOT submit
# with sbatch). It captures and recreates the software environment used by
# scripts 1-6, so the module's storage/memory habits are paired with an
# equally important habit: recording exactly which tool VERSIONS produced a
# given result.
#
# NEW IDEA on top of script 6: everything so far assumed the HPC's
# `module load` system was available. That works great on THIS cluster, but
# it is not portable — a labmate on a different cluster, a cloud VM, or their
# own laptop has no `module` command at all, and even on this same cluster,
# "SRA-Toolkit" might silently point to a newer default version next year.
# Conda/Mamba environments solve this: an environment.yml file is a portable,
# version-pinned recipe that recreates the *same* toolchain anywhere Conda or
# Mamba is installed, cluster or not.
#
# Usage (run once per fresh environment, or whenever the toolchain changes):
#   bash 7_Environment_Reproducibility.sh export   # snapshot current env -> environment.yml
#   bash 7_Environment_Reproducibility.sh create   # build a new env FROM environment.yml
#   bash 7_Environment_Reproducibility.sh verify   # confirm the active env matches environment.yml
###############################################################################

set -o pipefail

ENV_NAME="mibo8110-module0"
ENV_FILE="environment.yml"
ACTION="${1:-}"

# Sanity-check that some form of Conda/Mamba is actually available before
# doing anything else — the same "check your assumptions before acting"
# habit used for srr_list.txt in scripts 2-6, just applied to tooling
# instead of input data.
if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA_BIN="conda"
else
    echo "ERROR: neither 'mamba' nor 'conda' was found on PATH."
    echo "On most HPC systems this means you need to load a Conda/Mamba"
    echo "module first, e.g.: module load Miniconda3 (check 'module avail')."
    exit 1
fi
echo "Using '${CONDA_BIN}' as the environment manager."

write_environment_yml() {
    # A hand-written, version-pinned recipe covering every tool this module's
    # scripts rely on (1-6 plus the format-validation and capstone scripts
    # that follow). Pinning exact versions — not just package names — is the
    # entire point: "fastp" alone can mean a different tool behavior next
    # year, but "fastp=0.23.4" means the same behavior everywhere, forever.
    cat > "$ENV_FILE" <<'YAML'
name: mibo8110-module0
channels:
  - bioconda
  - conda-forge
dependencies:
  - sra-tools=3.0.3
  - fastqc=0.11.9
  - fastp=0.23.4
  - multiqc=1.28
  - samtools=1.16.1
  - bcftools=1.16
  - git=2.43
YAML
    echo "Wrote a version-pinned recipe to ${ENV_FILE}."
}

case "$ACTION" in
    export)
        # If the named environment already exists and is active-able, export
        # its ACTUAL resolved state (exact build strings, full dependency
        # tree) rather than the hand-written recipe above. This is the
        # "measure, don't guess" instinct from script 5/6 applied to
        # software versions: what you THINK you installed and what actually
        # got resolved can differ once transitive dependencies are involved.
        if "$CONDA_BIN" env list | grep -q "^${ENV_NAME} "; then
            "$CONDA_BIN" env export -n "$ENV_NAME" --no-builds > "$ENV_FILE"
            echo "Exported the ACTUAL resolved environment '${ENV_NAME}' to ${ENV_FILE}."
        else
            echo "No environment named '${ENV_NAME}' exists yet to export."
            echo "Writing a hand-authored starting recipe instead — build it with:"
            echo "  bash $0 create"
            write_environment_yml
        fi
        ;;
    create)
        if [ ! -s "$ENV_FILE" ]; then
            echo "${ENV_FILE} not found — writing the default recipe first."
            write_environment_yml
        fi
        echo "Creating environment '${ENV_NAME}' from ${ENV_FILE} (this can take a few minutes)..."
        "$CONDA_BIN" env create -n "$ENV_NAME" -f "$ENV_FILE"
        echo "Done. Activate it with: conda activate ${ENV_NAME}"
        ;;
    verify)
        if [ ! -s "$ENV_FILE" ]; then
            echo "ERROR: ${ENV_FILE} not found. Run '$0 export' or '$0 create' first."
            exit 1
        fi
        if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
            echo "WARNING: no Conda environment is currently active."
            echo "Run 'conda activate ${ENV_NAME}' first, then re-run verify."
            exit 1
        fi
        echo "Active environment: ${CONDA_DEFAULT_ENV}"
        if [ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]; then
            echo "WARNING: active environment name does not match '${ENV_NAME}' recorded in ${ENV_FILE}."
        fi
        # Spot-check a few key tool versions against what environment.yml
        # promises, the same way script 9 will spot-check FILE structure
        # against what a format promises.
        echo "Installed key tool versions in this environment:"
        "$CONDA_BIN" list -n "$CONDA_DEFAULT_ENV" 2>/dev/null | grep -E '^(sra-tools|fastqc|fastp|multiqc|samtools|bcftools|git)\b' || \
            echo "  (none of the expected tools were found — environment may be incomplete)"
        ;;
    *)
        echo "Usage: bash $0 {export|create|verify}"
        echo "  export  Snapshot the current/active Conda environment to ${ENV_FILE}"
        echo "  create  Build a new environment named '${ENV_NAME}' from ${ENV_FILE}"
        echo "  verify  Confirm the active environment matches what ${ENV_FILE} promises"
        exit 1
        ;;
esac

###############################################################################
# NUANCES & PITFALLS
#
#   - `module load` vs Conda are NOT the same kind of reproducibility.
#     `module load SAMtools/1.16.1-GCC-11.3.0` reproducibly gives you that
#     tool version ONLY on clusters that happen to host that exact module —
#     it means nothing on a laptop or a different center's cluster.
#     environment.yml is portable across any machine with Conda/Mamba
#     installed. Good practice on HPC is often BOTH: `module load` for
#     cluster-provided compilers/MPI, Conda for the bioinformatics toolchain
#     itself, recorded in environment.yml either way.
#
#   - Pinning package names without versions (e.g. just "samtools" instead
#     of "samtools=1.16.1") silently lets the resolved version drift every
#     time the environment is rebuilt — defeating the entire purpose of an
#     environment.yml. Always pin versions for anything that touches results.
#
#   - `conda env export` without `--no-builds` embeds OS/build-specific
#     strings (e.g. compiler ABI hashes) that can make the file fail to
#     resolve on a different OS/architecture. `--no-builds` trades a little
#     precision for a lot more portability — usually the right trade for a
#     teaching or collaboration environment.yml.
#
#   - This script's OWN environment.yml is exactly the kind of file script 8
#     (Git-based reproducibility) will commit to version control next —
#     an environment spec that is never committed alongside the code that
#     needs it is not actually reproducible, no matter how carefully it was
#     written.
###############################################################################
