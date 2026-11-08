#!/usr/bin/env bash
# Finite paired-read routing regression for BBCMS. Author: Sucrose.
# Run in a fresh result directory: bash dev/bloommultipasspairingtest.sh outdir=PATH [-ea|-da]
# Keep all inputs, output files, stderr and exit codes for baseline/fixed comparison.
set -euo pipefail
outdir=
assertions=-ea
for arg in "$@"; do
 case "$arg" in
  outdir=*) outdir=${arg#outdir=};;
  -ea|-da) assertions=$arg;;
  *) printf 'Unknown argument: %s\n' "$arg" >&2; exit 2;;
 esac
done
if [[ -z $outdir || -e $outdir ]]; then
 printf 'Supply outdir=PATH naming a directory that does not exist.\n' >&2
 exit 2
fi
tree=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p -- "$outdir"
outdir=$(cd -- "$outdir" && pwd)
# Four distinct 32-base reads. Mates have no reverse-complement overlap, so the
# first pass's default merge option cannot merge them. Quality 20 is unambiguous
# Phred+33 even though multipass does not forward qin/qout to the first pass.
printf '%s\n' '@pairA/1' 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' '+' '55555555555555555555555555555555' \
 '@pairB/1' 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT' '+' '55555555555555555555555555555555' > "$outdir/in1.fq"
printf '%s\n' '@pairA/2' 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' '+' '55555555555555555555555555555555' \
 '@pairB/2' 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG' '+' '55555555555555555555555555555555' > "$outdir/in2.fq"
common=("in=$outdir/in1.fq" "in2=$outdir/in2.fq" k=31 k0=31 bits=4 bits0=4
 cells=2000003 t=1 seed=0 ordered=t ecc0=f tossjunk0=f mincount0=0
 ecc=f merge=f ecco=f markerrors=f tossjunk=f mincount=0 minprob=0
 loglog=f pigz=f unpigz=f qin=33 qout=33 ow=f)
printf 'case\tprogram_exit\tfirst_mate_cmp\tsecond_mate_cmp\n' > "$outdir/checks.tsv"
failed=0
for name in one-pass-explicit two-pass-auto two-pass-explicit; do
 case "$name" in
  one-pass-explicit) mode=(passes=1 interleaved=f);;
  two-pass-auto) mode=(passes=2 interleaved=auto);;
  two-pass-explicit) mode=(passes=2 interleaved=f);;
 esac
 command=(bash "$tree/bbcms.sh" -Xmx256m -Xms32m "$assertions" "${common[@]}"
  "${mode[@]}" "out=$outdir/$name.1.fq" "out2=$outdir/$name.2.fq")
 printf '%q ' "${command[@]}" > "$outdir/$name.command"
 printf '\n' >> "$outdir/$name.command"
 code=0; first=0; second=0
 "${command[@]}" </dev/null > "$outdir/$name.stdout" 2> "$outdir/$name.stderr" || code=$?
 cmp "$outdir/in1.fq" "$outdir/$name.1.fq" > "$outdir/$name.first.diff" 2>&1 || first=$?
 cmp "$outdir/in2.fq" "$outdir/$name.2.fq" > "$outdir/$name.second.diff" 2>&1 || second=$?
 printf '%s\t%s\t%s\t%s\n' "$name" "$code" "$first" "$second" >> "$outdir/checks.tsv"
 if ((code!=0 || first!=0 || second!=0)); then failed=1; fi
done
# Additional controls retain the original three cases above verbatim.
# Paired files need not have names recognizable as adjacent interleaved mates.
sed 's,pairA/1,left.alpha,;s,pairB/1,left.beta,' "$outdir/in1.fq" > "$outdir/unusual1.fq"
sed 's,pairA/2,right.omega,;s,pairB/2,right.delta,' "$outdir/in2.fq" > "$outdir/unusual2.fq"
cat "$outdir/in1.fq" > "$outdir/unpaired.fq"
head -n 4 "$outdir/in2.fq" >> "$outdir/unpaired.fq"
# Remove only the two input arguments from the shared configuration.
settings=("${common[@]:2}")
for name in unusual-one-pass unusual-two-pass unusual-two-pass-alias unpaired-one-pass unpaired-two-pass-alias; do
 case "$name" in
  unusual-one-pass) mode=(passes=1 interleaved=f); paired_case=true;;
  unusual-two-pass) mode=(passes=2 interleaved=f); paired_case=true;;
  unusual-two-pass-alias) mode=(passes=2 int=f); paired_case=true;;
  unpaired-one-pass) mode=(passes=1 interleaved=f); paired_case=false;;
  unpaired-two-pass-alias) mode=(passes=2 int=f); paired_case=false;;
 esac
 if "$paired_case"; then
  input=("in=$outdir/unusual1.fq" "in2=$outdir/unusual2.fq")
  output=("out=$outdir/$name.1.fq" "out2=$outdir/$name.2.fq")
  expected=$outdir/unusual1.fq
 else
  input=("in=$outdir/unpaired.fq")
  output=("out=$outdir/$name.1.fq")
  expected=$outdir/unpaired.fq
 fi
 command=(bash "$tree/bbcms.sh" -Xmx256m -Xms32m "$assertions" "${input[@]}"
  "${settings[@]}" "${mode[@]}" "${output[@]}")
 printf '%q ' "${command[@]}" > "$outdir/$name.command"
 printf '\n' >> "$outdir/$name.command"
 code=0; first=0; second=-
 "${command[@]}" </dev/null > "$outdir/$name.stdout" 2> "$outdir/$name.stderr" || code=$?
 cmp "$expected" "$outdir/$name.1.fq" > "$outdir/$name.first.diff" 2>&1 || first=$?
 if "$paired_case"; then
  second=0
  cmp "$outdir/unusual2.fq" "$outdir/$name.2.fq" > "$outdir/$name.second.diff" 2>&1 || second=$?
  if ((second!=0)); then failed=1; fi
 fi
 printf '%s\t%s\t%s\t%s\n' "$name" "$code" "$first" "$second" >> "$outdir/checks.tsv"
 if ((code!=0 || first!=0)); then failed=1; fi
done
cat "$outdir/checks.tsv"
if ((failed)); then
 printf 'BLOOM_MULTIPASS_PAIRING_TEST_FAIL: inspect retained files in %s\n' "$outdir" >&2
 exit 1
fi
printf 'BLOOM_MULTIPASS_PAIRING_TEST_PASS: 8 cases, 14 exact output-file comparisons.\n'
