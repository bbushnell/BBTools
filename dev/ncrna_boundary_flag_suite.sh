#!/usr/bin/env bash
# Dori-only packaged-resource fixtures. Every case starts a fresh JVM.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
argument=${1:?outdir=/absolute/new/fixture/directory}
[[ "$argument" = outdir=/* ]]
out=${argument#outdir=}
test ! -e "$out"
mkdir -p "$out"
printf '>fixture\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$out/input.fa"
base=("in=$out/input.fa" t=1 taxonomy=f)
pass(){
  local name=$1
  local expectedp=$2
  local expectedt=$3
  shift 3
  bash "$root/dev/ncrna_boundary_flag_test.sh" "expectrnasep=$expectedp" "expecttmrna=$expectedt" "$@" "${base[@]}" > "$out/$name.stdout" 2> "$out/$name.stderr"
  grep -F 'BOUNDARY_FLAG_PASS ' "$out/$name.stdout"
}
fail(){
  local name=$1
  local message=$2
  shift 2
  if bash "$root/callgenes.sh" "${base[@]}" "$@" > "$out/$name.stdout" 2> "$out/$name.stderr";then
    printf 'Expected loud failure: %s\n' "$name" >&2
    exit 1
  fi
  grep -F "$message" "$out/$name.stderr"
}
pass defaults t t
pass explicit_off f f expectsmall=f expectlarge=f expect13=f expect85=f rnasepboundarynn=f tmrnaboundarynn=f srpsmallboundarynn=f srplargeboundarynn=f sixsrf00013boundarynn=f sixsrf01685boundarynn=f
pass rnasep_only t f rnasepboundarynn=t tmrnaboundarynn=f
pass tmrna_only f t tmrnaboundarynn=t rnasepboundarynn=f
pass both t t ncrnaboundarynet=f rnasepboundarynn=t tmrnaboundarynn=t
pass reverse_order t t rnasepboundarynn=t tmrnaboundarynn=t ncrnaboundarynet=f
pass srp_small t t expectsmall=t expectlarge=f srpsmallboundarynn=t srplargeboundarynn=f
pass srp_large t t expectsmall=f expectlarge=t srpsmallboundarynn=f srplargeboundarynn=t
pass srp_both t t expectsmall=t expectlarge=t srpsmallboundarynn=t srplargeboundarynn=t ncrnaboundarynet=f
pass ncrna_disabled f f expectsmall=f expectlarge=f expect13=f expect85=f expectmask=0 ncrna=f
pass tmrna_default_disabled t f expectmask=61 tmrna=f
fail rnasep_disabled 'rnasepboundarynn=t requires ncrna=t' ncrna=f rnasepboundarynn=t
fail tmrna_disabled 'tmrnaboundarynn=t requires ncrna=t and tmrna=t' tmrna=f tmrnaboundarynn=t
fail srpsmall_disabled 'srpsmallboundarynn=t requires ncrna=t' ncrna=f srpsmallboundarynn=t
fail srplarge_disabled 'srplargeboundarynn=t requires ncrna=t' ncrna=f srplargeboundarynn=t
fail global_missing_sixs 'Missing required sixs_rf00013 resource boundary:' ncrnaboundarynet=t "sixsrf00013boundarynet=$out/missing.bbnet"
printf 'BOUNDARY_FLAG_SUITE_PASS positive=11 negative=5\n'
