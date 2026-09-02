#!/usr/bin/env bash
###############################################################################
# 7_Environment_Reproducibility_updated.sh
#
# Purpose:
#   Create, verify, and document the Conda/Mamba software environment used by
#   this module.
#
# Where to keep and run it:
#   Keep this script and environment.yml in the persistent project repository
#   under HOME. The script may be launched from any working directory because
#   it resolves environment.yml relative to its own location.
#
# Do not submit this helper with sbatch. Run it interactively on the login node.
# Store large sequencing data and replaceable intermediate files in SCRATCH.
#
# Expected layout:
#   module-00-foundations-reproducibility/
#   `-- scripts/
#       |-- 7_Environment_Reproducibility_updated.sh
#       `-- environment.yml
#
# Usage:
#   bash 7_Environment_Reproducibility_updated.sh create
#   conda activate mibo8110-module0
#   bash 7_Environment_Reproducibility_updated.sh verify
#   bash 7_Environment_Reproducibility_updated.sh export
#   bash 7_Environment_Reproducibility_updated.sh versions
###############################################################################

set -Eeuo pipefail

readonly ENV_NAME="mibo8110-module0"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly ENV_FILE="${SCRIPT_DIR}/environment.yml"
readonly RESOLVED_ENV_FILE="${SCRIPT_DIR}/environment.resolved.yml"
readonly ACTION="${1:-}"

readonly EXPECTED_PACKAGES=(
  sra-tools
  fastqc
  fastp
  multiqc
  samtools
  bcftools
  git
)

log()  { printf '[INFO] %s\n' "$*"; }
warn() { printf '[WARNING] %s\n' "$*" >&2; }
die()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

log "Script directory: ${SCRIPT_DIR}"
log "Environment recipe: ${ENV_FILE}"

case "${SCRIPT_DIR}" in
  "${HOME}"|"${HOME}"/*)
    log "The script is stored under HOME. This is appropriate for persistent project files."
    ;;
  /scratch/*|/gpfs/scratch/*|/lustre/scratch/*|/tmp/*)
    warn "The script is under temporary or scratch storage. Keep the authoritative copy under HOME."
    ;;
esac

# Find an available environment manager.
if command -v mamba >/dev/null 2>&1; then
  CONDA_BIN="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BIN="conda"
else
  die "Neither mamba nor conda is available. Load the appropriate HPC module, then rerun the script. For example, inspect available modules with: module avail"
fi
readonly CONDA_BIN
log "Using ${CONDA_BIN}: $(command -v "${CONDA_BIN}")"

environment_exists() {
  "${CONDA_BIN}" env list 2>/dev/null |
    awk -v env_name="${ENV_NAME}" '
      $1 == env_name { found=1 }
      END { exit(found ? 0 : 1) }
    '
}

require_recipe() {
  [[ -s "${ENV_FILE}" ]] || die "Missing or empty recipe: ${ENV_FILE}. Place environment.yml beside this script."
}

show_versions() {
  local target_env="${1:-${ENV_NAME}}"

  environment_exists || die "Environment '${ENV_NAME}' does not exist. Run this script with 'create' first."

  log "Key package versions in '${target_env}':"
  "${CONDA_BIN}" list -n "${target_env}" 2>/dev/null |
    awk 'BEGIN {IGNORECASE=1}
      !/^#/ && $1 ~ /^(sra-tools|fastqc|fastp|multiqc|samtools|bcftools|git)$/ {
        printf "  %-12s %s\n", $1, $2
      }'
}

case "${ACTION}" in
  create)
    require_recipe

    if environment_exists; then
      log "Environment '${ENV_NAME}' already exists. Updating it from ${ENV_FILE}."
      "${CONDA_BIN}" env update \
        --name "${ENV_NAME}" \
        --file "${ENV_FILE}" \
        --prune
    else
      log "Creating environment '${ENV_NAME}' from ${ENV_FILE}."
      "${CONDA_BIN}" env create \
        --name "${ENV_NAME}" \
        --file "${ENV_FILE}"
    fi

    log "Activate the environment with: conda activate ${ENV_NAME}"
    ;;

  verify)
    require_recipe
    environment_exists || die "Environment '${ENV_NAME}' does not exist. Run this script with 'create' first."

    if [[ "${CONDA_DEFAULT_ENV:-}" != "${ENV_NAME}" ]]; then
      die "Environment '${ENV_NAME}' is not active. Run: conda activate ${ENV_NAME}"
    fi

    missing=0
    for package in "${EXPECTED_PACKAGES[@]}"; do
      if ! "${CONDA_BIN}" list -n "${ENV_NAME}" 2>/dev/null |
           awk -v package="${package}" '$1 == package {found=1} END {exit(found ? 0 : 1)}'; then
        warn "Missing expected package: ${package}"
        missing=1
      fi
    done

    show_versions "${ENV_NAME}"
    (( missing == 0 )) || die "Verification failed because one or more expected packages are missing."
    log "Verification completed successfully."
    ;;

  export)
    environment_exists || die "Environment '${ENV_NAME}' does not exist. Run this script with 'create' first."

    # Preserve the concise, curated environment.yml. Write the complete solved
    # package list to a separate reproducibility record instead.
    "${CONDA_BIN}" env export \
      --name "${ENV_NAME}" \
      --no-builds |
      sed '/^prefix:/d' > "${RESOLVED_ENV_FILE}"

    log "Exported the resolved environment to ${RESOLVED_ENV_FILE}."
    log "The curated recipe ${ENV_FILE} was not overwritten."
    ;;

  versions)
    show_versions "${ENV_NAME}"
    ;;

  *)
    cat <<USAGE
Usage: bash $(basename "$0") {create|verify|export|versions}

  create    Create or update '${ENV_NAME}' from environment.yml
  verify    Confirm that '${ENV_NAME}' is active and contains expected tools
  export    Write the solved package list to environment.resolved.yml
  versions  Display versions of the key tools
USAGE
    exit 1
    ;;
esac
