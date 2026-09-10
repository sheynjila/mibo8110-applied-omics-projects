# Module 3 — Microbial Genome Assembly: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- `spades.py` (tested with v3.15.5) and `flye` (tested with v2.9.6) —
  scripts 2 and 3 load these automatically via
  [`scripts/load_modules.sh`](scripts/load_modules.sh). Unlike other
  modules, this module's pinned versions were never verified against a
  real cluster's module tree (the original scripts assumed these tools
  were already on `PATH` with no `module load` at all) — run `module
  spider SPAdes`/`module spider Flye` yourself and set
  `SPADES_MODULE`/`FLYE_MODULE` if the pinned build isn't found.
- `minimap2` and `samtools` on `PATH` — used by script 6's read-depth check
  (also available via `load_modules.sh`'s `load_minimap2`/`load_samtools`;
  script 6 is plain Python and can't source `load_modules.sh` itself, so
  load these in your shell **before** running script 6). Script 6 degrades
  gracefully with a warning and skips the depth check if either is
  missing; the GC-content check still runs.
- Python 3 (no third-party packages required — scripts 4-7 use only the
  standard library, deliberately, so this module has no dependency chain
  beyond the assemblers themselves).
- At least ~9 GB of free RAM if you plan to run Flye **without** `--low-mem`
  on real bacterial genome-scale data (see the memory note below); with
  `--low-mem`, a few hundred MB is enough.
- Module 1 completed or reviewed — this module assumes the same
  "validate/interpret before you trust it" habit taught there, applied one
  level up the pipeline (assembled contigs instead of raw reads).

## A real, verified memory constraint (read before running Flye on a small machine)

Flye's `--nano-raw` mode defaults to k-mer size 17 for its k-mer counting
step, which allocates a fixed **~8.6 GB** flat array (`4^17/2` bytes)
**regardless of your genome's actual size** — a 200 kb toy genome and a 5 Mb
real bacterial genome both trigger the same fixed allocation. On any machine
with less than ~9 GB of free RAM, the default settings will crash with
`std::bad_alloc` inside `KmerCounter::count`. Use script 3's `--low-mem` flag
(k=15, ~537 MB) on such machines. On a real HPC node with the RAM to spare,
prefer the tuned default (k=17) for better repeat resolution — `--low-mem` is
a memory-constrained workaround, not a strict improvement.

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 1 | `bash 1_Choose_Assembly_Strategy.sh --reads1 reads_1.fastq [--reads2 reads_2.fastq]` | Console recommendation: `SHORT_READ`, `LONG_READ`, `HYBRID`, or `UNCERTAIN`, based on sampled read lengths |
| 2a | (if SHORT_READ/HYBRID) `bash 2_Run_SPAdes_Short_Read_Assembly.sh reads_1.fastq reads_2.fastq spades_out --threads 8 --mem 32` | `spades_out/contigs.fasta`, `spades_out/assembly_graph_with_scaffolds.gfa` |
| 2b | (if LONG_READ/HYBRID) `bash 3_Run_Flye_Long_Read_Assembly.sh reads.fastq flye_out --threads 8 [--low-mem]` | `flye_out/assembly.fasta`, `flye_out/assembly_graph.gfa`, `flye_out/assembly_info.txt` |
| 3 | `python3 4_Interpret_Assembly_Graph.py <out_dir>/assembly_graph*.gfa` | Console report: branch-point count (repeat signature) and circularity signature |
| 4 | `python3 5_Compute_Assembly_Statistics.py <out_dir>/<contigs\|assembly>.fasta --label "<assembler>"` | Console report: N50/L50/N90/L90, per-contig length+GC%, mandatory N50 caveat |
| 5 | `python3 6_Assess_Completeness_And_Contamination.py <out_dir>/<contigs\|assembly>.fasta --reads1 reads_1.fastq [--reads2 reads_2.fastq]` | Console report: per-contig GC and depth outlier flags, combined verdict |
| 6 | `python3 7_Assembly_Quality_Report_And_Fitness_For_Purpose.py <out_dir>/<contigs\|assembly>.fasta <out_dir>/assembly_graph*.gfa --reads1 reads_1.fastq [--reads2 reads_2.fastq] --assembler "<name>" --out ASSEMBLY_REPORT.md` | `ASSEMBLY_REPORT.md` written, sections 1-3 auto-filled, sections 4-5 blank |
| 7 | Fill in Section 4 (interpretation) and Section 5 (fitness-for-purpose conclusion) of `ASSEMBLY_REPORT.md` by hand | Completed portfolio artifact |

Steps 3-6 must run in order for a single assembler's output — each is
independent enough to run standalone, but script 7 calls scripts 4-6
internally as subprocesses to build the report, so all three must be able to
run cleanly against the same `<out_dir>` first.

If script 1 recommends `HYBRID`, run **both** 2a and 2b against the same
organism's short- and long-read files respectively, then run steps 3-6
**twice** (once per assembler output) to compare — this module does not
include a dedicated read-combining hybrid assembler mode; it teaches
strategy selection and independent evaluation of each assembler's result,
which is itself a valid, common first pass before attempting a true hybrid
assembly with a tool like Unicycler.

## Expected outputs, by step

- **Step 1:** stdout only (no new files) — strategy recommendation.
- **Step 2a:** `spades_out/contigs.fasta`, `spades_out/scaffolds.fasta`,
  `spades_out/assembly_graph_with_scaffolds.gfa`, `spades_out/spades.log`.
- **Step 2b:** `flye_out/assembly.fasta`, `flye_out/assembly_graph.gfa`,
  `flye_out/assembly_info.txt`, `flye_out/flye.log`.
- **Step 3:** stdout only — topology report.
- **Step 4:** stdout only — statistics report.
- **Step 5:** stdout only — completeness/contamination report.
- **Step 6:** `ASSEMBLY_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Script 1 prints `UNCERTAIN` | Sampled read lengths don't clearly fall into either the short-read or long-read range (e.g. a very small or unusual test file) | Inspect the printed length distribution; if genuinely mixed/ambiguous, treat as HYBRID and run both assembler paths for comparison |
| Flye crashes with `std::bad_alloc` (often inside a stack trace mentioning `KmerCounter` or `vertex_index`) | Machine has less than ~9 GB free RAM and Flye is using its default k=17 k-mer counting, which allocates a fixed ~8.6 GB array regardless of genome size | Re-run with script 3's `--low-mem` flag (k=15, ~537 MB) — this is a documented, verified workaround, not a bug in your data |
| SPAdes runs but reports far more contigs than expected for a small/simple genome | Read coverage too low for the genome's actual complexity (very short reads, low read count, or a low-error but under-covered simulation) | Increase read count/coverage; SPAdes generally needs at least ~20-30x depth across the genome to assemble cleanly, more if the genome has extended repeats |
| Script 4 reports 0 branch points and 0 circular contigs on a graph you expected to show one or the other | Some assemblers report circularity only in a separate summary file (Flye's `assembly_info.txt` `circ.` column), not always as a literal self-loop `L` line in every GFA variant | Cross-check `assembly_info.txt`/`contigs.paths` alongside script 4's GFA-only view rather than treating script 4 as the sole source of truth for circularity |
| Script 4 reports a branch point on an assembly that fully collapsed to one contig in the final FASTA | Expected and informative, not a bug — SPAdes/Flye can resolve a repeat at the scaffolding stage (using paired-end or long-read spanning evidence) even though the underlying assembly graph still shows the ambiguous branch point where the repeat itself sits. The final contig is one piece; the graph still honestly records that a repeat was present and had to be resolved | Report this exactly as it appears: contiguous final assembly, with a graph-level note that a repeat was detected and spanned — this is stronger evidence than a graph with no repeat signal at all |
| Script 6 prints "minimap2/samtools not found on PATH — skipping the read-depth outlier check" | Those tools are not installed | Install them, or accept GC-only contamination screening (weaker but still informative) |
| Script 6 flags your single, correct contig as an "outlier" and does *not* flag known-contaminant contigs | You are running an out-of-date copy of the script using an unweighted per-contig median/mean; a genome fragmented into many small pieces plus one large correct contig can invert the result under an unweighted average | Use the current version of script 6, which weights both the GC reference and the depth median by contig **length**, not contig count |
| Script 6 reports "No contigs flagged by either proxy check" | Either genuinely no proxy-detectable contamination, or contamination too subtle for GC/depth alone to catch | Report as "no proxy evidence of contamination" (as the script itself states), not "confirmed clean" — run real BUSCO/CheckM if a stronger completeness/contamination claim is needed |
| Script 7's `ASSEMBLY_REPORT.md` sections 4-5 are blank | Expected — these are intentionally left for the analyst to complete by hand | Fill them in based on the actual topology (script 4), statistics (script 5), and contamination findings (script 6), and a stated intended use for the assembly |

## Completion checklist (portfolio artifact)

The curriculum's Module 3 assessment asks for: *"Assembly FASTA, environment
and commands, assembly statistics, graph or coverage evidence, quality
interpretation, and a clear fitness-for-purpose conclusion."* Confirm each is
present before considering the module complete:

- [ ] **Assembly FASTA** — the actual `contigs.fasta`/`assembly.fasta` from
      step 2a/2b is retained alongside the report, not just summarized.
- [ ] **Environment and commands** — the exact assembler, version, and
      command-line flags used (including whether `--low-mem`/`--hybrid` was
      set) are recorded — copy them from your step 2a/2b invocation.
- [ ] **Assembly statistics** — Section 2 of `ASSEMBLY_REPORT.md` reflects a
      real script 5 run against your actual assembly, including the N50
      caveat text, not a placeholder.
- [ ] **Graph/coverage evidence** — Section 1 reflects a real script 4 run
      (branch points and circularity), and Section 3 reflects a real script
      6 run (GC/depth outlier flags), both against your actual data.
- [ ] **Quality interpretation** — Section 4 is filled in by hand,
      referencing the actual topology, statistics, and contamination
      findings above — not left as the auto-generated prompt text.
- [ ] **Fitness-for-purpose conclusion** — Section 5 is filled in by hand,
      stating explicitly what the assembly will and will not be used for and
      whether the evidence supports that specific use.
