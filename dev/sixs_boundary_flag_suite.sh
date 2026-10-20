#!/usr/bin/env bash
# CLI-to-loaded6S configuration checks; borrowed RNaseP resources test plumbing only.
# No genomic accuracy claim is made with these deliberately borrowed resources.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
argument=${1:?outdir=/absolute/new/fixture/directory}
[[ "$argument" = outdir=/* ]]
out=${argument#outdir=}
test ! -e "$out"
mkdir -p "$out"
printf '>fixture\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$out/input.fa"
base=("in=$out/input.fa" t=1 taxonomy=f)
resources13=("sixsrf00013boundarynet=$root/resources/callgenes_RNaseP/boundary.bbnet" "sixsrf00013boundarystarttable=$root/resources/callgenes_RNaseP/tables/start.tsv" "sixsrf00013boundarystoptable=$root/resources/callgenes_RNaseP/tables/stop.tsv")
resources85=("sixsrf01685boundarynet=$root/resources/callgenes_RNaseP/boundary.bbnet" "sixsrf01685boundarystarttable=$root/resources/callgenes_RNaseP/tables/start.tsv" "sixsrf01685boundarystoptable=$root/resources/callgenes_RNaseP/tables/stop.tsv")
positive=0
negative=0
pass(){
  local name=$1
  shift
  bash "$root/dev/ncrna_boundary_flag_test.sh" expectrnasep=t expecttmrna=t "$@" "${base[@]}" > "$out/$name.stdout" 2> "$out/$name.stderr"
  grep -F 'SIXS_BOUNDARY_FLAG_PASS ' "$out/$name.stdout"
  ((positive+=1))
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
  ((negative+=1))
}
pass defaults expect13=t expect85=t
pass packaged expect13=t expect85=t sixsrf00013boundarynn=t sixsrf01685boundarynn=t
pass rf13 expect13=t expect85=f sixsrf00013boundarynn=t sixsrf01685boundarynn=f "${resources13[@]}"
pass rf85 expect13=f expect85=t sixsrf00013boundarynn=f sixsrf01685boundarynn=t "${resources85[@]}"
pass both expect13=t expect85=t "${resources13[@]}" "${resources85[@]}" sixsrf00013boundarynn=t sixsrf01685boundarynn=t
pass global expect13=t expect85=t ncrnaboundarynet=t "${resources13[@]}" "${resources85[@]}"
pass optout13 expect13=f expect85=t sixsrf00013boundarynn=f ncrnaboundarynet=t "${resources85[@]}"
pass optout85 expect13=t expect85=f ncrnaboundarynet=t sixsrf01685boundarynn=f "${resources13[@]}"
pass euk_gate expect13=t expect85=t euk=t
pass ranges expect13=t expect85=t expect13start=-2,-1,0,1 expect13stop=0,1,2 expect85start=-1,0 expect85stop=0 \
  "${resources13[@]}" "${resources85[@]}" sixsrf00013boundarystartrange=-2,1 sixsrf00013boundarystoprange=0,2 \
  sixsrf01685boundarystartrange=-1,0 sixsrf01685boundarystoprange=0,0 sixsrf00013boundarynn=t sixsrf01685boundarynn=t
for family in rf00013 rf01685;do
  prefix=sixs$family
  fail "${family}_parent" "require ncrna=t and sixs=t" ncrna=f "${prefix}boundarynn=t"
  fail "${family}_sixs" "require ncrna=t and sixs=t" sixs=f "${prefix}boundarynn=t"
  fail "${family}_resources_off" "require effective boundary refinement enabled" "${prefix}boundarynn=f" "${prefix}boundarynet=not-opened"
  fail "${family}_range_off" "require effective boundary refinement enabled" "${prefix}boundarynn=f" "${prefix}boundarystartrange=-1,1"
  fail "${family}_path_parent" "require ncrna=t and sixs=t" ncrnaboundarynet=t sixs=f "${prefix}boundarynet=not-opened"
  fail "${family}_range_parent" "require ncrna=t and sixs=t" ncrnaboundarynet=t sixs=f "${prefix}boundarystoprange=-1,1"
  fail "${family}_bad_range" "Boundary range must include zero" "${prefix}boundarynn=t" "${prefix}boundarystartrange=1,2"
  fail "${family}_missing" "Missing required sixs_${family} resource boundary:" "${prefix}boundarynn=t" "${prefix}boundarynet=$out/missing.bbnet"
done
fail swapped_start 'boundary start table has type STOP' sixsrf00013boundarynn=t "${resources13[@]}" "sixsrf00013boundarystarttable=$root/resources/callgenes_RNaseP/tables/stop.tsv"
fail swapped_stop 'boundary stop table has type START' sixsrf01685boundarynn=t "${resources85[@]}" "sixsrf01685boundarystoptable=$root/resources/callgenes_RNaseP/tables/start.tsv"
test "$positive" = 10
test "$negative" = 18
printf 'SIXS_BOUNDARY_FLAG_SUITE_PASS positive=%s negative=%s scope=configuration_only\n' "$positive" "$negative"
