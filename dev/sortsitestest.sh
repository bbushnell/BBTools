#!/usr/bin/env bash
# Literal two-site CLI checks for ordering and caller-selected temporary filenames.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
bbtools=$(cd "$(dirname "$SCRIPT")/.." && pwd)
out=; phase=fixed
for arg in "$@"; do
	case "$arg" in
	outdir=*) out=${arg#outdir=};;
	bbtools=*) bbtools=${arg#bbtools=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	*) echo 'Use outdir=new-directory [bbtools=absolute-root] [phase=baseline|fixed]' >&2; exit 1;;
	esac
done
[[ $bbtools = /* && -f $bbtools/javasetup.sh && -n $out && ! -e $out && -d ${out%/*} ]]
mkdir -- "$out"
out=$(cd "$out" && pwd)
DIR=$bbtools
source "$bbtools/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment
first='1,0,10,13,4,2,0,00,0,100'
second='1,0,40,43,4,1,0,00,0,100'
printf '%s\t%s\n' "$first" "$second" > "$out/input.txt"
printf '%s\n%s\n' "$second" "$first" > "$out/id.expected"
printf '%s\t%s\n' "$first" "$second" > "$out/position.expected"
printf 'case\texpected\tobserved\tresult\n'
for case in default-id custom-id custom-position; do
	mkdir "$out/$case"
	mode=id; key=0
	if [[ $case = custom-position ]]; then mode=position; key=4294967296; fi
	flags=()
	if [[ $case != default-id ]]; then flags=("tempname=$out/$case/custom_#.txt"); fi
	command=(java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$bbtools/current"
		pacbio.SortSites "$out/input.txt" "$out/$case/output.txt" "mode=$mode" deletefiles=f "${flags[@]}")
	printf '%q ' "${command[@]}" > "$out/$case/command.txt"
	printf '\n' >> "$out/$case/command.txt"
	(cd "$out/$case" && "${command[@]}") > "$out/$case/stdout.txt" 2> "$out/$case/stderr.txt"
	cmp "$out/$mode.expected" "$out/$case/output.txt"
	printf '%s-output\texact_sorted_records\texact\tPASS\n' "$case"
	legacy=$out/$case/SortSitesByIDTempFile_b-1_$key.txt.gz
	custom=$out/$case/custom_b-1_$key.txt
	if [[ $case = default-id ]]; then
		test -f "$legacy"
		printf '%s-path\tdefault_bucket\tdefault_bucket\tPASS\n' "$case"
	elif [[ $phase = baseline ]]; then
		test -f "$legacy" && test ! -e "$custom"
		printf '%s-path\tignored_custom_path\tdefault_bucket\tCAPTURED\n' "$case"
	else
		test -f "$custom" && test ! -e "$legacy"
		printf '%s-path\tcustom_bucket\tcustom_bucket\tPASS\n' "$case"
	fi
done
