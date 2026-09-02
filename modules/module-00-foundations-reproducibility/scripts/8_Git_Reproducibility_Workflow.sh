#!/bin/bash
###############################################################################
# 8_Git_Reproducibility_Workflow.sh
#
# INTERACTIVE helper — run this directly on the login node (do NOT submit
# with sbatch). It sets up (or verifies) a Git repository for a project
# directory and embeds the exact commit hash into a run's log output, so
# that any set of results can always be traced back to the precise version
# of the code that produced them.
#
# NEW IDEA on top of script 7: script 7 pinned WHICH TOOLS and versions were
# used. This script pins WHICH CODE was used. A tool version alone doesn't
# tell you whether the script itself was edited between two runs that
# produced different numbers — only a commit hash does that. Together,
# environment.yml (script 7) + a Git commit hash (this script) answer the
# two questions every reproducibility check eventually asks: "what software,
# and what code, produced this result?"
#
# Usage:
#   bash 8_Git_Reproducibility_Workflow.sh init      # one-time repo setup
#   bash 8_Git_Reproducibility_Workflow.sh snapshot  # commit current state
#   bash 8_Git_Reproducibility_Workflow.sh stamp     # print a run-log header with the current commit
###############################################################################

set -o pipefail

ACTION="${1:-}"

if ! command -v git >/dev/null 2>&1; then
    echo "ERROR: 'git' was not found on PATH."
    echo "Load it first, e.g.: module load git   (or install it via script 7's Conda environment)."
    exit 1
fi

write_gitignore() {
    # Mirrors the repository-wide .gitignore categories already established
    # for this project: never commit sequencing data, Slurm logs, or other
    # multi-gigabyte pipeline OUTPUTS — only commit CODE, specs, and small
    # text artifacts. Committing a 3GB FASTQ file is the single most common
    # way students accidentally make a Git repository unusable.
    cat > .gitignore <<'IGNORE'
# --- HPC / Slurm job artifacts ---
slurm-*.out
logs/
*.log
mem_usage_log.txt

# --- Sequencing data and large intermediate/output files ---
*.fastq
*.fastq.gz
*.fasta
*.fa
*.fa.gz
*.fna
*.fna.gz
*.gff
*.gff.gz
*.gtf
*.sam
*.bam
*.bai
*.vcf
*.vcf.gz
*.sra
srr_list.txt

# --- Conda ---
.conda/

# --- OS / editor noise ---
.DS_Store
*.swp
IGNORE
    echo "Wrote .gitignore (matching the repository's existing data-safety categories)."
}

case "$ACTION" in
    init)
        if [ -d .git ]; then
            echo "A Git repository already exists here ($(pwd)). Nothing to initialize."
        else
            git init
            git config user.name  "${GIT_AUTHOR_NAME:-$(whoami)}"
            git config user.email "${GIT_AUTHOR_EMAIL:-$(whoami)@example.edu}"
            echo "Initialized a new Git repository in $(pwd)."
        fi
        [ -f .gitignore ] || write_gitignore
        echo "Next: add your scripts and environment.yml, then run '$0 snapshot'."
        ;;
    snapshot)
        if [ ! -d .git ]; then
            echo "ERROR: no Git repository here yet. Run '$0 init' first."
            exit 1
        fi
        [ -f .gitignore ] || write_gitignore
        # Stage everything NOT excluded by .gitignore. `-A` respects the
        # .gitignore written above, so raw sequencing data, Slurm logs, and
        # other pipeline OUTPUTS are never staged in the first place.
        #
        # (An earlier version of this script tried to stage only specific
        # extensions with `git add -- '*.sh' '*.md' ...`, but git treats a
        # pathspec with zero matches as a hard error and aborts the ENTIRE
        # add -- so if a project had no *.R files yet, nothing at all got
        # staged, silently. `-A .` avoids that trap.)
        git add -A .
        # Second, explicit line of defense: warn (don't just trust the
        # .gitignore) if anything large or clearly data-shaped slipped
        # through -- e.g. a --force-added file, or an extension the
        # .gitignore doesn't yet know about.
        BIG_STAGED=$(git diff --cached --name-only -z | xargs -0 -I{} du -m "{}" 2>/dev/null | awk '$1 > 20 {print}')
        if [ -n "$BIG_STAGED" ]; then
            echo "WARNING: the following staged file(s) are over 20MB -- confirm they"
            echo "         are really meant to be committed, not pipeline output:"
            echo "$BIG_STAGED"
        fi
        if git diff --cached --quiet; then
            echo "Nothing new to commit — working tree already matches the last snapshot."
        else
            git commit -m "Snapshot: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
            echo "Committed a new snapshot: $(git rev-parse --short HEAD)"
        fi
        ;;
    stamp)
        if [ ! -d .git ]; then
            echo "ERROR: no Git repository here yet. Run '$0 init' first."
            exit 1
        fi
        # This is the piece meant to be embedded directly into another
        # script's log output (see script 10, which calls this in `stamp`
        # mode and redirects it into its own run log) — a one-line,
        # unambiguous answer to "what exact code produced this run?"
        if ! git diff --quiet || ! git diff --cached --quiet; then
            echo "WARNING: uncommitted changes are present — this run's code state was"
            echo "         NOT fully captured by a commit. Run '$0 snapshot' before"
            echo "         submitting jobs whose results you plan to keep."
        fi
        echo "git_commit=$(git rev-parse HEAD 2>/dev/null || echo 'NO_COMMITS_YET')"
        echo "git_branch=$(git branch --show-current 2>/dev/null || echo 'unknown')"
        echo "git_dirty=$(git diff --quiet && git diff --cached --quiet && echo false || echo true)"
        ;;
    *)
        echo "Usage: bash $0 {init|snapshot|stamp}"
        echo "  init      One-time: create the repository and a data-safe .gitignore"
        echo "  snapshot  Commit the current code + environment.yml state"
        echo "  stamp     Print commit/branch/dirty-state lines for embedding in a run log"
        exit 1
        ;;
esac

###############################################################################
# NUANCES & PITFALLS
#
#   - A commit hash only proves what it claims if the working tree was
#     CLEAN at run time. `stamp` above deliberately warns (rather than
#     silently stamping) when uncommitted changes exist — a run log that
#     says "commit abc123" while the actual code had un-committed edits is
#     actively misleading, worse than no version information at all.
#
#   - Git tracks CODE well and DATA poorly. The .gitignore written here
#     blocks the same large-file categories as scripts elsewhere in this
#     project for a reason: a single accidentally-committed FASTQ file can
#     bloat a repository permanently, even after the file is later deleted,
#     because Git retains every version in its history by default.
#
#   - `git commit` records a snapshot of tracked files, not proof that the
#     snapshot was ever actually RUN. Pair every meaningful result with
#     both a commit hash (this script) AND a log showing the command that
#     was executed (scripts 1-6's own stdout/mem_usage_log.txt) — hash
#     alone answers "what code existed," not "what code ran."
#
#   - Committing on the login node is appropriate (it is a lightweight,
#     interactive, metadata-only operation) — unlike the compute-heavy work
#     in scripts 1-6 and 10, it does not belong in an sbatch job.
###############################################################################
