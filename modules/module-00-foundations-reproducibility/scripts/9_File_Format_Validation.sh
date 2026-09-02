#!/bin/bash
###############################################################################
# 9_File_Format_Validation.sh
#
# INTERACTIVE helper — run this directly on the login node (do NOT submit
# with sbatch). It scans a directory for the file formats this curriculum
# relies on most (FASTA, FASTQ, SAM/BAM, BED, GFF/GTF, VCF) and runs a
# lightweight structural sanity check on each one it finds, printing a
# short PASS/WARN/FAIL line per file.
#
# NEW IDEA on top of scripts 7-8: those two scripts pinned WHICH TOOLS and
# WHICH CODE produced a result. This script checks whether the DATA those
# tools consumed and produced actually has the shape it claims to. A tool
# given a truncated FASTQ file or a header-less VCF will often still run —
# and produce silently wrong output — rather than crash outright. These
# checks are cheap enough to run before AND after every pipeline step.
#
# Usage:
#   bash 9_File_Format_Validation.sh <directory>
#
# Example:
#   bash 9_File_Format_Validation.sh capstone/module1/raw_reads
###############################################################################

set -o pipefail

TARGET_DIR="${1:-.}"

if [ ! -d "$TARGET_DIR" ]; then
    echo "ERROR: '${TARGET_DIR}' is not a directory."
    echo "Usage: bash $0 <directory>"
    exit 1
fi

PASS=0
WARN=0
FAIL=0

report() {
    # $1=status $2=file $3=message
    case "$1" in
        PASS) PASS=$((PASS + 1)) ;;
        WARN) WARN=$((WARN + 1)) ;;
        FAIL) FAIL=$((FAIL + 1)) ;;
    esac
    printf "%-4s  %-55s  %s\n" "$1" "$2" "$3"
}

check_fastq() {
    local f="$1"
    local reader="cat"
    case "$f" in *.gz) reader="zcat" ;; esac
    local n_lines
    n_lines=$($reader "$f" 2>/dev/null | wc -l)
    if [ "$n_lines" -eq 0 ]; then
        report FAIL "$f" "empty file"
        return
    fi
    # A valid FASTQ file's line count must be a multiple of 4: header,
    # sequence, '+' separator, quality — one full record every 4 lines.
    if [ "$((n_lines % 4))" -ne 0 ]; then
        report FAIL "$f" "line count (${n_lines}) is not a multiple of 4 -- truncated file?"
        return
    fi
    local first_char
    first_char=$($reader "$f" 2>/dev/null | head -1 | cut -c1)
    if [ "$first_char" != "@" ]; then
        report FAIL "$f" "first record does not start with '@' -- not a FASTQ header"
        return
    fi
    report PASS "$f" "$((n_lines / 4)) records, header format OK"
}

check_fastq_pair() {
    # Paired-end FASTQ files must have IDENTICAL read counts, or every
    # downstream paired-end tool (fastp, aligners) will either error out or,
    # worse, silently mis-pair reads.
    local r1="$1" r2="$2"
    local reader1="cat" reader2="cat"
    case "$r1" in *.gz) reader1="zcat" ;; esac
    case "$r2" in *.gz) reader2="zcat" ;; esac
    local n1 n2
    n1=$($reader1 "$r1" 2>/dev/null | wc -l); n1=$((n1 / 4))
    n2=$($reader2 "$r2" 2>/dev/null | wc -l); n2=$((n2 / 4))
    if [ "$n1" -ne "$n2" ]; then
        report FAIL "$(basename "$r1") + $(basename "$r2")" "read counts differ (${n1} vs ${n2}) -- pairing broken"
    else
        report PASS "$(basename "$r1") + $(basename "$r2")" "${n1} reads in both mates -- pairing OK"
    fi
}

check_fasta() {
    local f="$1"
    local reader="cat"
    case "$f" in *.gz) reader="zcat" ;; esac
    local n_headers
    n_headers=$($reader "$f" 2>/dev/null | grep -c '^>')
    if [ "$n_headers" -eq 0 ]; then
        report FAIL "$f" "no '>' header lines found -- not a valid FASTA"
    else
        report PASS "$f" "${n_headers} sequence header(s) found"
    fi
}

check_gff_gtf() {
    local f="$1"
    local reader="cat"
    case "$f" in *.gz) reader="zcat" ;; esac
    # A GFF/GTF file has 9 tab-separated columns per feature line (comments
    # starting with '#' are excluded). Tallying feature TYPES (column 3)
    # is a fast way to sanity-check that the annotation contains what a
    # pipeline expects (e.g. "gene", "exon", "CDS") before trusting it.
    local n_features
    n_features=$($reader "$f" 2>/dev/null | grep -v '^#' | awk -F'\t' 'NF==9' | wc -l)
    if [ "$n_features" -eq 0 ]; then
        report FAIL "$f" "no 9-column feature lines found -- malformed or empty"
        return
    fi
    report PASS "$f" "${n_features} feature lines; top types: $($reader "$f" 2>/dev/null | grep -v '^#' | awk -F'\t' 'NF==9 {print $3}' | sort | uniq -c | sort -rn | head -3 | awk '{printf "%s(%s) ", $2, $1}')"
}

check_bed() {
    local f="$1"
    # BED requires at least 3 tab/space-separated columns (chrom, start,
    # end) with start < end. Checking this catches the classic "0-based
    # vs 1-based off-by-one" mistake before it silently corrupts an
    # intersect/overlap analysis downstream.
    local bad
    bad=$(awk 'NF<3 {c++} NF>=3 && $2>=$3 {c++} END{print c+0}' "$f")
    local total
    total=$(wc -l < "$f")
    if [ "$bad" -gt 0 ]; then
        report FAIL "$f" "${bad}/${total} lines have <3 columns or start>=end"
    else
        report PASS "$f" "${total} intervals, all start<end"
    fi
}

check_vcf() {
    local f="$1"
    local reader="cat"
    case "$f" in *.gz) reader="zcat" ;; esac
    local header_line
    header_line=$($reader "$f" 2>/dev/null | head -1)
    if [[ "$header_line" != "##fileformat=VCF"* ]]; then
        report FAIL "$f" "missing/invalid '##fileformat=VCF...' first line"
        return
    fi
    local col_header
    col_header=$($reader "$f" 2>/dev/null | grep -m1 '^#CHROM')
    if [ -z "$col_header" ]; then
        report FAIL "$f" "no '#CHROM' column-header line found"
        return
    fi
    local n_variants
    n_variants=$($reader "$f" 2>/dev/null | grep -vc '^#')
    report PASS "$f" "valid header, ${n_variants} variant record(s)"
}

check_sam_bam() {
    local f="$1"
    if ! command -v samtools >/dev/null 2>&1; then
        report WARN "$f" "samtools not on PATH -- skipped (load it via script 7's environment)"
        return
    fi
    # `samtools quickcheck` is purpose-built for exactly this: a fast
    # integrity check (valid header, EOF block present for BAM) without
    # reading every alignment record.
    if samtools quickcheck "$f" 2>/dev/null; then
        local n_reads
        n_reads=$(samtools view -c "$f" 2>/dev/null)
        report PASS "$f" "samtools quickcheck OK, ${n_reads} alignment record(s)"
    else
        report FAIL "$f" "samtools quickcheck failed -- truncated or corrupt"
    fi
}

echo "Scanning '${TARGET_DIR}' for known bioinformatics file formats..."
echo ""

# FASTA
while IFS= read -r -d '' f; do check_fasta "$f"; done < <(find "$TARGET_DIR" -maxdepth 2 \( -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fa.gz' -o -iname '*.fasta.gz' -o -iname '*.fna' -o -iname '*.fna.gz' \) -print0)

# FASTQ (individual check, then pair up files ending in _1/_2 or _R1/_R2)
mapfile -d '' FASTQS < <(find "$TARGET_DIR" -maxdepth 2 \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) -print0)
for f in "${FASTQS[@]}"; do check_fastq "$f"; done
for r1 in "${FASTQS[@]}"; do
    case "$r1" in
        *_1.fastq*|*_R1.fastq*|*_1.fq*|*_R1.fq*)
            r2="${r1/_1./_2.}"; r2="${r2/_R1./_R2.}"
            [ -f "$r2" ] && check_fastq_pair "$r1" "$r2"
            ;;
    esac
done

# SAM/BAM
while IFS= read -r -d '' f; do check_sam_bam "$f"; done < <(find "$TARGET_DIR" -maxdepth 2 \( -iname '*.bam' -o -iname '*.sam' \) -print0)

# BED
while IFS= read -r -d '' f; do check_bed "$f"; done < <(find "$TARGET_DIR" -maxdepth 2 -iname '*.bed' -print0)

# GFF/GTF
while IFS= read -r -d '' f; do check_gff_gtf "$f"; done < <(find "$TARGET_DIR" -maxdepth 2 \( -iname '*.gff' -o -iname '*.gff3' -o -iname '*.gtf' -o -iname '*.gff.gz' -o -iname '*.gtf.gz' \) -print0)

# VCF
while IFS= read -r -d '' f; do check_vcf "$f"; done < <(find "$TARGET_DIR" -maxdepth 2 \( -iname '*.vcf' -o -iname '*.vcf.gz' \) -print0)

echo ""
echo "Summary: ${PASS} passed, ${WARN} warnings, ${FAIL} failed."
if [ "$FAIL" -gt 0 ]; then
    exit 1
fi

###############################################################################
# NUANCES & PITFALLS
#
#   - These checks are STRUCTURAL, not BIOLOGICAL. A FASTQ file can pass
#     every check here and still contain garbage sequence, and a VCF can
#     have a perfect header and still contain wrong genotype calls. Format
#     validation catches "the file is broken," not "the file is correct" —
#     it is a floor, not a ceiling, exactly like FastQC in script 1 tells
#     you about read quality but not about whether the sample itself was
#     mislabeled.
#
#   - `samtools quickcheck` deliberately does NOT decompress and verify
#     every single alignment record (that would defeat the "quick" in its
#     name). It catches truncated files and missing EOF markers reliably,
#     but a BAM with a few corrupted records in the middle can still pass.
#     For anything downstream that must be bulletproof, follow up with
#     `samtools view` piped through a full read, or `samtools stats`.
#
#   - The FASTQ record-count-multiple-of-4 check catches TRUNCATION but not
#     every corruption mode — a file with the exact wrong number of lines
#     inserted at a 4-line boundary would still "pass" this check while
#     being scrambled. It is a cheap, high-value first filter, not an
#     exhaustive validator.
#
#   - Run this BEFORE trusting a downloaded/received file, and again AFTER
#     any step that rewrites it (trimming, alignment, variant calling) —
#     the same "measure before and after" habit that scripts 1-6 apply to
#     read counts and resource usage, just applied to file structure.
###############################################################################
