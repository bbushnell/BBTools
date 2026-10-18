#!/usr/bin/env bash
# Raiden: ordered GFF IDs must continue across inputs, including no-call batches.
set -euo pipefail
export LC_ALL=C
root=$(cd "$(dirname "$0")/.." && pwd)
[[ "${1:-}" = outdir=* ]]
out=${1#outdir=}
test ! -e "$out"
mkdir -p "$out"
trap 'printf "%s\n" "$?" > "$out/exit_code.txt"' EXIT
gene=NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCTGGTGGCCATAGCGAGGTGGAAACACCCGTTCCCATCCCGAACACGGAAGTTAAGCCCCTCAGCGCCGATGGTACTGTGGGGTCACCCCGTGGGAGAGTAGGTCGCCGCCAGGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
printf '>a_gene\n%s\n>a_blank\nNNNNNNNNNNNNNNNNNNNN\n' "$gene" > "$out/a.fa"
printf '>b_blank\nNNNNNNNNNNNNNNNNNNNN\n' > "$out/b.fa"
printf '>c_gene1\n%s\n>c_gene2\n%s\n' "$gene" "$gene" > "$out/c.fa"
flags=(5sonly=t taxonomy=f readbufferlength=1 t=1 -Xmx1g ow=f)
for name in a b c;do
  bash "$root/callgenes.sh" "in=$out/$name.fa" "out=$out/$name.gff" "${flags[@]}" ordered=t > "$out/$name.stdout" 2> "$out/$name.stderr"
done
cat "$out/a.gff" > "$out/expected.gff"
sed '/^#/d' "$out/b.gff" "$out/c.gff" >> "$out/expected.gff"
test "$(grep -vc '^#' "$out/expected.gff")" = 3
cut -f1 "$out/expected.gff" | grep -v '^#' > "$out/names.txt"
printf 'a_gene\nc_gene1\nc_gene2\n' | cmp - "$out/names.txt"
inputs=$out/a.fa,$out/b.fa,$out/c.fa
bash "$root/callgenes.sh" "in=$inputs" "out=$out/ordered.gff" "${flags[@]}" ordered=t t=2 > "$out/ordered.stdout" 2> "$out/ordered.stderr"
cmp "$out/expected.gff" "$out/ordered.gff"
bash "$root/callgenes.sh" "in=$inputs" "out=$out/unordered.gff" "${flags[@]}" ordered=f > "$out/unordered.stdout" 2> "$out/unordered.stderr"
sort "$out/ordered.gff" > "$out/ordered.sorted.gff"
sort "$out/unordered.gff" > "$out/unordered.sorted.gff"
cmp "$out/ordered.sorted.gff" "$out/unordered.sorted.gff"
printf 'PASS ordered three-input GFF, multiple lists, no-call batches, exact single-input concatenation and unordered content parity\n'
