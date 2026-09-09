# Master Runbook

This is the **program-level operational** document: environment setup,
how manual vs. orchestrated execution both work, and the cross-cutting
troubleshooting that applies to more than one module. For a specific
module's exact execution order and expected outputs, use that module's own
`RUNBOOK.md` — this document is the shared context those all build on.

## 1. Prerequisites (every module)

- A Slurm-managed HPC cluster (`sbatch`, `squeue`, `sacct`, `module` all on
  `PATH`) — every pipeline script targets Slurm. See
  [Computing environment guidance](README.md#computing-environment-guidance-per-curriculum-4)
  in the root README for smaller-scale alternatives (WSL, Galaxy) while
  learning concepts.
- Git, for cloning this repo and (if you're an instructor extending it)
  committing changes.
- Familiarity with Module 0's content — every later module assumes the
  shell, environment, and reproducibility habits taught there.

## 2. Getting the repo onto the cluster

```bash
git clone https://github.com/<owner>/mibo8110-applied-omics-projects.git
cd mibo8110-applied-omics-projects
```

Scripts assume they're run from *inside* their own module's `scripts/`
directory, or submitted with a relative path from the repo root (both
work — see §4 below for why).

## 3. Running a module manually

1. Read that module's `README.md` (the "why" — guiding question, learning
   outcomes, portfolio artifact) and then its `RUNBOOK.md` (the "how" —
   exact commands, expected outputs, troubleshooting).
2. `cd` into that module's `scripts/` directory (or note its path).
3. Stage the inputs the RUNBOOK asks for — almost always an `srr_list.txt`
   with one SRA accession per line, in the working directory the script
   will create under `/scratch/$(whoami)/...`.
4. Submit the script: `sbatch <script>.sh` for anything with an `#SBATCH`
   header, or `bash <script>.sh [args]` for a plain validation/analysis
   script (check the script's own header comment — modules 1, 4, and 5 mix
   both kinds).
5. Follow the module's RUNBOOK execution-order table step by step.

## 4. Environment-module loading — how `load_modules.sh` works

Every pipeline script that needs an external tool (`fastqc`, `bwa`,
`kraken2`, ...) sources a `load_modules.sh` living alongside it in the same
`scripts/` directory, instead of calling `module load <exact-build>`
directly. Three things to know:

- **Try pinned → discover → refuse.** [`lib/module_loader.sh`](lib/module_loader.sh)
  tries the exact build this course was authored against first. If that's
  missing on your cluster, it searches `module avail <tool-name>`: exactly
  one match gets loaded (with a loud `WARNING` — check your run log for
  these), and zero or multiple matches is a hard, explicit `ERROR` rather
  than a silent guess.
- **Override per tool** with `<LABEL>_MODULE=<exact-name>` (label is the
  lowercase name in that module's `load_modules.sh`, e.g. `fastqc` →
  `FASTQC_MODULE`). Find the exact name yourself with `module spider
  <tool>` if the automatic discovery refuses.
- **Why sourcing, not executing:** `module load` only changes the
  environment of the *current* shell. Every `load_modules.sh` refuses to
  run standalone (`bash load_modules.sh` prints an error) — it must be
  `source`d from the pipeline script that needs it.

### The SLURM script-staging gotcha (why path resolution looks the way it does)

Every retrofitted script resolves its own directory like this:

```bash
SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-NN-name/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"
```

not just `dirname "${BASH_SOURCE[0]}"`. Under real `sbatch`, a submitted
script can be staged/copied into a spool location, where `BASH_SOURCE`
resolves to a path that does **not** contain the script's sibling files
(like `load_modules.sh`) — a bare `dirname "$0"` lookup then fails to find
them. `$SLURM_SUBMIT_DIR` (the directory `sbatch` was invoked from) is
reliable in that situation; `$PIPELINE_SCRIPT_DIR` is an explicit override
the [orchestration copies](pipeline_orchestration/README.md) set so they
can point back at the original module's `load_modules.sh` without
duplicating it.

**Practical implication:** submit with `sbatch` from the repo root (so
`$SLURM_SUBMIT_DIR` + the known relative path resolves correctly), or from
inside the script's own `scripts/` directory (so the plain `dirname`
fallback resolves correctly). Either works; a `cd` to some unrelated
directory before submitting does not.

## 5. Manual vs. orchestrated execution

See [`pipeline_orchestration/README.md`](pipeline_orchestration/README.md)
for the full explanation. Short version: run the scripts under
`modules/*/scripts/` by hand while learning; use the `RUN_TAG`-parameterized
copies under `pipeline_orchestration/scripts/` when you need many
concurrent/repeated runs that must not collide on one working directory.
Both call the same underlying tools with the same flags — the only
difference is `WORKDIR` isolation and where `load_modules.sh` is sourced
from.

## 6. Cross-cutting troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ERROR: no module found for <tool>` | Neither the pinned build nor anything matching `module avail <tool>` exists on this cluster | Run `module spider <tool>` yourself, then re-run with `<LABEL>_MODULE=<exact-name>` |
| `ERROR: <tool> is ambiguous on this cluster` | Multiple modules matched the loose search term | Pick the right one from the printed list and set `<LABEL>_MODULE=<exact-name>` |
| `load_modules.sh must be sourced, not executed` | Ran `bash load_modules.sh` or `./load_modules.sh` directly | `source` it from a pipeline script instead — `module load` has no effect once the sourcing shell exits |
| MultiQC fails or picks up the wrong Python | Another module (R, QIIME2, a pinned tool's own bundled Python) won the ambiguous `python3` resolution before MultiQC ran | Every module that calls MultiQC loads `load_python` immediately before `load_multiqc`, in that order — if you added a new MultiQC call, keep that order |
| A pipeline script can't find `load_modules.sh` / other sibling files under `sbatch` | Classic SLURM script-staging path issue — see §4 above | Submit from the repo root or from inside that `scripts/` directory, not an unrelated `cwd` |
| Long `git commit -m "..."` silently truncates special characters (`$`, `` ` ``, `#`) | Bash/shell interpolation inside a quoted `-m` argument | Write the message to a file and commit with `git commit -F <file>`; verify with `git log -1 --format='%B'` |
| A multi-sample loop reports one sample "FAILED" but the whole job kept going | Working as intended — every multi-sample script wraps each sample in its own subshell precisely so one failure doesn't abort the run | Check that sample's own log/output for the real cause; the job's overall exit code is not the signal to look at |

## 7. For instructors extending this repo

- Follow the two-tier README/RUNBOOK pattern (root README.md/root
  RUNBOOK.md are the program view of a page mirrored, per module, and
  every existing module folder). New content belongs in a module's own
  folder, linked from that module's README — see the root README's
  "README architecture" section.
- New environment-module dependencies get one `load_<tool>()` line added
  to that module's `scripts/load_modules.sh` (or `lib/module_loader.sh`
  directly, for the shared engine), never a bare `module load` in the
  pipeline script itself.
- A new orchestratable master pipeline gets a copy under
  `pipeline_orchestration/scripts/<module>/`, changed only by the
  `RUN_TAG`/`PIPELINE_SCRIPT_DIR` pattern shown in that directory's
  existing copies — never edited to add orchestration logic in place.
