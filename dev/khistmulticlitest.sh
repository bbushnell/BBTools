#!/bin/bash
# Small CLI checks. Set KHIST_TEST_THREADS=4 on a suitable host for thread parity.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
work=$(mktemp -d /tmp/khistmulti-cli.XXXXXX)
printf 'KHISTMULTI_CLI_DIR %s\n' "$work"
bash "$root/dev/khistmultitest.sh" "$work"
common=(k=32,64,96,124 buckets=2048 seed=19 minprob=0 raw=t qin=33 -Xmx256m)
bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" int=t t=1 out="$work/paired.tsv"
bash "$root/khistmulti.sh" "${common[@]}" in="$work/r1.fq" in2="$work/r2.fq" t=1 out="$work/split.tsv"
cmp "$work/paired.tsv" "$work/split.tsv"
bash "$root/khistmulti.sh" "${common[@]}" in=stdin.fq int=t t=1 < "$work/paired.fq" > "$work/stdin.tsv"
cmp "$work/paired.tsv" "$work/stdin.tsv"
bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" int=t t="${KHIST_TEST_THREADS:-1}" out="$work/threaded.tsv"
cmp "$work/paired.tsv" "$work/threaded.tsv"
bash "$root/khistmulti.sh" -Xmx256m config="$root/testdata/khistmulti_test.config" \
  in="$work/paired.fq" out="$work/config.tsv" outstream="$work/config.diagnostics"
cmp "$work/paired.tsv" "$work/config.tsv"
grep -q 'K values:' "$work/config.diagnostics"
bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" int=t t=1 reads=0 out="$work/empty.tsv"
if grep -q '^[0-9]' "$work/empty.tsv"; then exit 1; fi
bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" int=t t=1 k=33,95,127,300 out="$work/arbitrary.tsv"
grep -q 'Count_K33' "$work/arbitrary.tsv"
for index in 1 2 3 4; do
  ks=(32 64 96 124)
  k=${ks[index-1]}
  bash "$root/khistmulti.sh" "${common[@]}" k="$k" in="$work/paired.fq" int=t t=1 out="$work/k$k.tsv"
  awk -v column="$((index*2))" '!/^#/ {print $1, $column, $(column+1)}' "$work/paired.tsv" > "$work/multi$k.txt"
  awk '!/^#/ {print $1, $2, $3}' "$work/k$k.tsv" > "$work/single$k.txt"
  diff -u "$work/multi$k.txt" "$work/single$k.txt"
done
if bash "$root/khistmulti.sh" "${common[@]}" in="$root/testdata/khistmulti_bad.fq" \
  int=f t=1 out="$work/malformed.tsv"; then
  printf 'Malformed FASTQ unexpectedly succeeded.\n' >&2
  exit 1
fi
test ! -e "$work/malformed.tsv"
for bad in k=0 k=32,32 histmax=65536 buckets=0 minprob=NaN unrecognized=t; do
  if bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" t=1 "$bad" out="$work/bad.tsv"; then
    printf 'Invalid argument unexpectedly succeeded: %s\n' "$bad" >&2
    exit 1
  fi
done
if bash "$root/khistmulti.sh" "${common[@]}" in="$work/paired.fq" t=1 out="$work/paired.tsv" overwrite=f; then
  printf 'Overwrite protection failed.\n' >&2
  exit 1
fi
printf 'KHISTMULTI_CLI_PASS\n'
