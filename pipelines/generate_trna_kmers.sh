#!/bin/bash
# generate_trna_kmers.sh — Reproduce the tRNA long-kmer screening set.
#
# Provenance (Aug 16, 2026):
#   Input: 2,515,952 tRNAs from two sources:
#     1. bacteria4+archaea4 corpus (577,687 tRNAs, 11,100 genomes)
#        extracted via annotateanticodon.sh + cutgff.sh types=tRNA minlen=40
#     2. NCBI RefSeq PGAP (gtdb232) corpus (1,938,265 tRNAs, ~36,600 genomes)
#        extracted via cutgff.sh types=tRNA minlen=40
#   Method: iterative greedy set-cover with kmercountexact.sh + bbduk.sh.
#     Each round finds the most common k-mers in the still-unmatched tRNAs,
#     then removes matched tRNAs. Descent cutoffs chosen to balance kmer count
#     vs coverage (minimize kmers, maximize detection).
#   Result: 2,868 13-mers covering 99.85% of all input tRNAs.
#     FP rate on random 75bp genomic shreds: <0.46% (includes real tRNA shreds;
#     true FP is lower).
#   Scope: prokaryotic tRNAs only. Chloroplast, mitochondrial, and fungal tRNAs
#     are NOT represented — add those corpora and re-run when available.
#
# Usage: generate_trna_kmers.sh in=<all_trnas.fa> out=<kmers.fa> [k=13]
#   in=  concatenated FASTA of all tRNA sequences (any source/extraction)
#   out= output FASTA of selected k-mers (the screening set)
#   k=   kmer length (default 13)
#
# The script is deterministic: same input → same output.

set -euo pipefail

# Parse args
IN="" OUT="" K=13
for arg in "$@"; do
  case "$arg" in
    in=*) IN="${arg#in=}" ;;
    out=*) OUT="${arg#out=}" ;;
    k=*) K="${arg#k=}" ;;
  esac
done
if [ -z "$IN" ] || [ -z "$OUT" ]; then
  echo "Usage: generate_trna_kmers.sh in=<all_trnas.fa> out=<kmers.fa> [k=13]" >&2
  exit 1
fi

TOTAL=$(grep -c '>' "$IN")
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

CUTOFFS="16000 8000 4000 1000 500 200 100 50 25 20 15"
CUMULATIVE=""
UNMATCHED="$IN"
CUM_KMERS=0

echo "Input: $TOTAL tRNAs, k=$K" >&2
echo "Round  Cutoff  NewKmers  CumKmers  Unmatched  Coverage%" >&2

round=0
for cutoff in $CUTOFFS; do
  round=$((round+1))

  kmercountexact.sh in="$UNMATCHED" k=$K \
    out="$TMPDIR/kmers_r${round}.fa" mincount=$cutoff \
    2>"$TMPDIR/r${round}.err"

  new_kmers=$(grep -c '>' "$TMPDIR/kmers_r${round}.fa" 2>/dev/null || echo 0)

  if [ "$new_kmers" -eq 0 ]; then
    printf "  %2d    %5d    %6d    %8d    %8d    (no new kmers)\n" \
      "$round" "$cutoff" 0 "$CUM_KMERS" "$(grep -c '>' "$UNMATCHED")" >&2
    continue
  fi

  if [ -z "$CUMULATIVE" ]; then
    cp "$TMPDIR/kmers_r${round}.fa" "$TMPDIR/kmers_cum.fa"
  else
    cat "$CUMULATIVE" "$TMPDIR/kmers_r${round}.fa" > "$TMPDIR/kmers_cum.fa"
  fi
  CUMULATIVE="$TMPDIR/kmers_cum.fa"
  CUM_KMERS=$(grep -c '>' "$CUMULATIVE")

  bbduk.sh in="$IN" ref="$CUMULATIVE" k=$K mm=f \
    out="$TMPDIR/unmatched_r${round}.fa" \
    2>"$TMPDIR/r${round}_bbduk.err"

  UNMATCHED="$TMPDIR/unmatched_r${round}.fa"
  remaining=$(grep -c '>' "$UNMATCHED")
  coverage=$(echo "scale=3; ($TOTAL - $remaining) * 100 / $TOTAL" | bc)

  printf "  %2d    %5d    %6d    %8d    %8d    %s%%\n" \
    "$round" "$cutoff" "$new_kmers" "$CUM_KMERS" "$remaining" "$coverage" >&2
done

cp "$CUMULATIVE" "$OUT"
remaining=$(grep -c '>' "$UNMATCHED")
echo "" >&2
echo "Wrote $CUM_KMERS kmers to $OUT" >&2
echo "Coverage: $(echo "scale=3; ($TOTAL - $remaining) * 100 / $TOTAL" | bc)% of $TOTAL tRNAs" >&2
echo "Unmatched: $remaining" >&2
