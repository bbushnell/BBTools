#!/usr/bin/env bash
# Geometry plus fresh-JVM CLI-to-family-to-scavenger binding; no genomic accuracy claim.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
argument=${1:?outdir=/absolute/new/directory}
[[ "$argument" = outdir=/* ]]
out=${argument#outdir=}
test ! -e "$out"
mkdir -p "$out"
printf '>fixture\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$out/input.fa"
base=("in=$out/input.fa" t=1 taxonomy=f)
pass(){
  local name=$1
  local pad=$2
  shift 2
  bash "$root/dev/ncrna_collapsed_window_padding_test.sh" "expectpad=$pad" "${base[@]}" "$@" > "$out/$name.stdout" 2> "$out/$name.stderr"
  grep -Fx "NCRNA_COLLAPSED_PADDING_CLI_PASS pad=$pad gene_min=80 other_families_unchanged=true" "$out/$name.stdout"
  grep -Fx 'NCRNA_COLLAPSED_PADDING_TEST_PASS default_exact median168 retained_exact separate_clusters physical_edges invalid_controls' "$out/$name.stdout"
}
fail(){
  local name=$1
  local expected=$2
  shift 2
  if bash "$root/callgenes.sh" "${base[@]}" "$@" > "$out/$name.stdout" 2> "$out/$name.stderr";then
    printf 'Expected failure: %s\n' "$name" >&2;exit 1
  fi
  grep -F "$expected" "$out/$name.stderr"
}
pass defaults 75
pass explicit_zero 0 srpsmallcollapsedpad=0
pass padding75 75 srpsmallcollapsedpad=75
fail negative 'srpsmallcollapsedpad must be >=0' srpsmallcollapsedpad=-1
fail too_short 'Collapsed-window padding must be disabled or span the unchanged gene-length floor' srpsmallcollapsedpad=10
fail disabled 'srpsmallcollapsedpad requires the SRP-small caller to be enabled' ncrna=f srpsmallcollapsedpad=75
fail wrong_family 'srpsmallcollapsedpad requires the SRP-small caller to be enabled' ncrnafamily=rnasep srpsmallcollapsedpad=75
printf 'NCRNA_COLLAPSED_PADDING_SUITE_PASS positive=3 negative=4 production_default=75\n'
