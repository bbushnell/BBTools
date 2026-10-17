#!/usr/bin/env bash
# Literal selector/extractor fixtures, including sequential reservations.
set -euo pipefail
runtime=${1:?tested runtime}; root=${2:?fresh fixture directory}; baseline=${3:-}
[[ ! -e $root ]]; mkdir -p "$root"
printf '0\ta\t2\n' > "$root/original.tsv"
printf 'active_index\tfamily_id\trep_id\tx\tparent_family_id\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\n0\t0\ta\tx\tNA\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\n' > "$root/roster.tsv"
printf 'a\ta\na\ta1\nb\tb\nb\tb1\nb\tb2\nc\tc\nc\tc1\nc\tc2\nd\td\ne\te\n' > "$root/clusters.tsv"
printf 'active_index\tfamily_id\trep_id\tmember_count\n1\t1\tb\t3\n' > "$root/prior.tsv"
printf 'active_index\tfamily_id\trep_id\tmember_count\n2\t2\tc\t3\n3\t3\td\t1\n' > "$root/expected.tsv"
run(){
 local rt=$1 output=$2;shift 2
 bash "$rt/dev/hbmfamilyexpand.sh" mode=select in="$root/clusters.tsv" original="$root/original.tsv" \
  roster="$root/roster.tsv" count=2 expectedfamilies=5 expectedgenes=10 out="$output" "$@" -Xmx256m -eoom
}
run "$runtime" "$root/actual" reserved="$root/prior.tsv" > "$root/select.out" 2> "$root/select.err"
cmp "$root/actual.selected.tsv" "$root/expected.tsv"
run "$runtime" "$root/next" reserved="$root/prior.tsv" reserved="$root/actual.selected.tsv" count=1 > "$root/next.out" 2> "$root/next.err"
printf 'active_index\tfamily_id\trep_id\tmember_count\n4\t4\te\t1\n' > "$root/expected_next.tsv"
cmp "$root/next.selected.tsv" "$root/expected_next.tsv"
run "$runtime" "$root/plain" > "$root/plain.out" 2> "$root/plain.err"
if [[ -n $baseline ]]; then
 run "$baseline" "$root/baseline" > "$root/baseline.out" 2> "$root/baseline.err"
 for suffix in selected.tsv excluded.tsv summary.tsv; do cmp "$root/plain.$suffix" "$root/baseline.$suffix"; done
fi
for kind in duplicate id_gap active_gap wrong_count absent reversed; do
 case $kind in
 duplicate) args=(reserved="$root/prior.tsv" reserved="$root/prior.tsv");;
 reversed) args=(reserved="$root/actual.selected.tsv" reserved="$root/prior.tsv");;
 *)
  case $kind in id_gap) row=$'1\t9\tb\t3';; active_gap)row=$'4\t1\tb\t3';; wrong_count)row=$'1\t1\tb\t2';; absent)row=$'1\t1\tz\t3';; esac
  head -n1 "$root/prior.tsv" > "$root/$kind.tsv"; printf '%s\n' "$row" >> "$root/$kind.tsv"
  args=(reserved="$root/$kind.tsv");;
 esac
 if run "$runtime" "$root/reject_$kind" "${args[@]}" > "$root/$kind.out" 2> "$root/$kind.err"; then
  echo "Invalid $kind selection succeeded" >&2;exit 1
 fi
 [[ ! -e $root/reject_$kind.selected.tsv ]]
done
# Split-child ancestry excludes the original a cluster even when representative differs.
sed 's/^0\t0\ta\tx\tNA/0\t5\tchild\tx\t0/' "$root/roster.tsv" > "$root/split.tsv"
run "$runtime" "$root/split" roster="$root/split.tsv" > "$root/split.out" 2> "$root/split.err"
printf 'active_index\tfamily_id\trep_id\tmember_count\n1\t6\tb\t3\n2\t7\tc\t3\n' > "$root/expected_split.tsv"
cmp "$root/split.selected.tsv" "$root/expected_split.tsv"
printf 'BBGAP_REFERENCE_SELECTION_PASS\n' > "$root/PASS"
cat "$root/PASS"
