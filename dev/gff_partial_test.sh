#!/bin/bash
# Actual reader, RNA trainer and all three CutGff implementations on synthetic input.
# Requires full gff/prok package compilation; optional BBTOOLS_TEST_CLASSES overlays native classes.
set -eo pipefail
script_path="${BASH_SOURCE[0]}"
while [[ -L "$script_path" ]]; do
	script_dir="$(cd -- "$(dirname -- "$script_path")" && pwd)"
	script_path="$(readlink -- "$script_path")"
	[[ "$script_path" = /* ]] || script_path="$script_dir/$script_path"
done
DIR="$(cd -- "$(dirname -- "$script_path")/.." && pwd)/"
CP="${DIR}current/"
if [[ -n "${BBTOOLS_TEST_CLASSES:-}" ]]; then
	[[ "$BBTOOLS_TEST_CLASSES" = /* && -d "$BBTOOLS_TEST_CLASSES" ]] || { echo 'BBTOOLS_TEST_CLASSES must name an existing absolute directory' >&2; exit 1; }
	CP="$BBTOOLS_TEST_CLASSES:$CP"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=32m --mode=fixed "$@"
setEnvironment
set -u
fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/gff-partial.XXXXXXXX")
printf 'FIXTURE_DIR=%s\n' "$fixture_dir"
run_java(){ java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" "$@"; }
run_java gff.GffPartialTest prepare "$fixture_dir"
for cutter in CutGff CutGff2 CutGff_ST; do
	for mode in filtered all banned required; do
		extra=(banpartial=t)
		case "$mode" in
			all) extra=(banpartial=f);;
			banned) extra+=(banattributes=keepout);;
			required) extra+=(attributes=wanted);;
		esac
		output="$fixture_dir/$cutter.$mode.fna"
		run_java "gff.$cutter" "in=$fixture_dir/input.fna" "gff=$fixture_dir/input.gff" \
			"out=$output" type=rRNA t=1 ow=t "${extra[@]}" > "$fixture_dir/$cutter.$mode.log" 2>&1 || {
			cat "$fixture_dir/$cutter.$mode.log" >&2
			exit 1
		}
		run_java gff.GffPartialTest verify "$output" "$mode"
	done
done
echo 'PASS GFF_PARTIAL_SUITE'
