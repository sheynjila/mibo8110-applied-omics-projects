#!/bin/bash
###############################################################################
# lib/module_loader.sh
#
# Shared HPC environment-module loading helper for this curriculum.
#
# WHY THIS EXISTS: every module's scripts were originally written with bare
# `module purge; module load <exact-pinned-version>` calls. That is fine on
# the cluster the course was authored on, but a real HPC site's module tree
# almost never matches another site's build-string-for-build-string (e.g.
# `FastQC/0.11.9-Java-11` vs `FastQC/0.12.1-Java11`), so a bare `module load`
# on a different cluster fails the whole job with no useful next step for a
# student. `load_module_capability()` tries the pinned build first (what
# this course was tested against), and only if that's absent falls back to
# discovering the tool by name -- loading it automatically if there is
# exactly one candidate (with a loud warning), and refusing to guess if
# there are zero or several, because silently picking the wrong version is
# worse than stopping with a clear message.
#
# USAGE (from a module's own scripts/load_modules.sh -- see any module's
# load_modules.sh for a concrete example):
#
#   source ".../lib/module_loader.sh"
#   load_fastqc() { load_module_capability fastqc "FastQC/0.11.9-Java-11" "FastQC"; }
#
# Then, in a pipeline script:
#
#   module purge
#   load_fastqc || exit 1
#
# HOW TO OVERRIDE ON YOUR OWN CLUSTER: export <LABEL>_MODULE=<exact-name>
# before running, e.g. `FASTQC_MODULE=FastQC/0.12.3-GCCcore-13.2.0 sbatch ...`.
# This is checked before the pinned build and before discovery, so it always
# wins.
#
# THIS FILE MUST BE SOURCED, NOT EXECUTED -- `module load` only changes the
# environment of the current shell; running this as a subprocess would load
# modules into a shell that immediately exits, doing nothing for the caller.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: lib/module_loader.sh must be sourced, not executed (module load does not survive a subshell)." >&2
    echo "  Fix: source \"\$(dirname \"\${BASH_SOURCE[0]}\")/load_modules.sh\"" >&2
    exit 1
fi

# resolve_script_dir [caller-relative-path]
#
# Returns the directory a pipeline script is running from. Resolution order:
#   1. $PIPELINE_SCRIPT_DIR, if the caller (or an orchestration wrapper) has
#      already exported it -- always wins, use this when you know better.
#   2. $SLURM_SUBMIT_DIR joined with the given caller-relative path, if
#      running under `sbatch` -- SLURM_SUBMIT_DIR is the directory `sbatch`
#      was invoked from and is reliable even when SLURM stages/copies the
#      submitted script, unlike BASH_SOURCE/$0 (which can resolve to a spool
#      path that does not contain this script's sibling files, e.g.
#      load_modules.sh). See docs/curriculum -- this bit us for real once.
#   3. dirname of BASH_SOURCE of the calling script -- correct for direct,
#      interactive (`bash script.sh`) execution.
resolve_script_dir() {
    local caller_rel_dir="${1:-}"
    if [[ -n "${PIPELINE_SCRIPT_DIR:-}" ]]; then
        echo "${PIPELINE_SCRIPT_DIR}"
    elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -n "${caller_rel_dir}" ]]; then
        echo "${SLURM_SUBMIT_DIR}/${caller_rel_dir}"
    else
        (cd "$(dirname "${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}")" && pwd)
    fi
}

# load_module_capability <label> <pinned-module-spec> <module-avail-search-term>
#
# <label>  short lowercase identifier, e.g. "fastqc" -- also the override
#          env var name, upper-cased with "_MODULE" appended (FASTQC_MODULE).
# <pinned-module-spec>  the exact module this course was authored/tested
#          against, e.g. "FastQC/0.11.9-Java-11".
# <module-avail-search-term>  a loose term to search with `module avail` when
#          the pinned build is missing, e.g. "FastQC".
load_module_capability() {
    local label="$1" pinned="$2" search="$3"
    local override_var
    override_var="$(echo "${label}" | tr '[:lower:]' '[:upper:]')_MODULE"
    local override="${!override_var:-}"

    if [[ -n "${override}" ]]; then
        if module load "${override}" 2>/dev/null; then
            echo "[load_modules] ${label}: loaded override module '${override}' (\$${override_var})"
            return 0
        fi
        echo "[load_modules] ERROR: \$${override_var}='${override}' failed to load. Check 'module spider ${override}'." >&2
        return 1
    fi

    if module load "${pinned}" 2>/dev/null; then
        echo "[load_modules] ${label}: loaded pinned module '${pinned}'"
        return 0
    fi

    echo "[load_modules] WARNING: pinned module '${pinned}' not found for ${label} -- searching 'module avail ${search}'..." >&2
    local avail n
    avail="$(module -t avail "${search}" 2>&1 | grep -v '^-' | grep -v '^\s*$' || true)"
    n="$(echo "${avail}" | grep -c . || true)"

    if [[ "${n}" -eq 1 ]]; then
        local found="${avail}"
        echo "[load_modules] WARNING: using discovered module '${found}' instead of pinned '${pinned}' for ${label} -- results may differ from the pinned version; note this in your report." >&2
        module load "${found}"
        return $?
    elif [[ "${n}" -eq 0 ]]; then
        echo "[load_modules] ERROR: no module found for ${label}. Pinned '${pinned}' is missing and 'module avail ${search}' found nothing on this cluster." >&2
        echo "  Fix: find the right module name yourself (e.g. 'module spider ${search}') and re-run with ${override_var}=<exact-name>." >&2
        return 1
    else
        echo "[load_modules] ERROR: ${label} is ambiguous on this cluster -- 'module avail ${search}' found ${n} candidates:" >&2
        echo "${avail}" | sed 's/^/    /' >&2
        echo "  Refusing to guess which one. Re-run with ${override_var}=<exact-name-from-the-list-above>." >&2
        return 1
    fi
}
