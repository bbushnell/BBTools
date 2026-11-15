#!/usr/bin/env bash
# Actual CLI pair: exact normal output and fatal missing-sequence worker error.
# Baseline uses one well so producer+terminal cannot fill the four-well queue.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
out=; phase=fixed; bbtools=$DIR; program=reformat
for arg in "$@"; do
	case "$arg" in
	outdir=*) out=${arg#outdir=};;
	bbtools=*) bbtools=${arg#bbtools=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	program=reformat|program=finder) program=${arg#program=};;
	*) echo 'Use outdir=new-child-of-existing-directory [phase=baseline|fixed] [program=reformat|finder] [bbtools=absolute-root]' >&2; exit 1;;
	esac
done
[[ $bbtools = /* && -f "$bbtools/reformatpb.sh" ]] || { echo 'bbtools must be an absolute root containing reformatpb.sh' >&2; exit 1; }
[[ -n $out && ! -e $out && -d ${out%/*} ]] || { echo 'Output must be a new child of an existing directory' >&2; exit 1; }
mkdir -- "$out"
out=$(cd "$out" && pwd)
wells=16; [[ $phase != baseline ]] || wells=1
printf '@HD\tVN:1.6\tSO:unknown\n' > "$out/normal.sam"
printf '@HD\tVN:1.6\tSO:unknown\n' > "$out/invalid.sam"
: > "$out/normal.expected.fq"
# Contiguous wells with three four-base reads each; first invalid interior has SEQ/QUAL=*.
for ((well=42; well<42+wells; well++)); do
	for ((part=0; part<3; part++)); do
		start=$((part*4)); stop=$((start+4))
		name=m/$well/${start}_${stop}
		printf '%s\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n' "$name" >> "$out/normal.sam"
		if [[ $well = 42 && $part = 1 ]]; then
			printf '%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n' "$name" >> "$out/invalid.sam"
		else
			printf '%s\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n' "$name" >> "$out/invalid.sam"
		fi
		printf '@%s\nACGT\n+\nIIII\n' "$name" >> "$out/normal.expected.fq"
	done
done
launcher=reformatpb.sh; class=ReformatPacBio; diagnostic='internal median length 0'
flags=(trim=f)
if [[ $program = finder ]]; then
	launcher=icecreamfinder.sh; class=IceCreamFinder; diagnostic=NullPointerException
	flags=(trim=t adapter=null alignrc=f alignadapter=f)
fi
normal=(bash "$bbtools/$launcher" "in=$out/normal.sam" "out=$out/normal.fq"
	threads=1 readbufferlength=1 readbuffers=2 "${flags[@]}" trimpolya=f minlen=0
	flaglongreads=f jni=f simd=f pigz=f unpigz=f qin=33 qout=33 ow=f -Xmx128m)
invalid=(bash "$bbtools/$launcher" "in=$out/invalid.sam" out=null
	threads=1 readbufferlength=1 readbuffers=2 "${flags[@]}" trimpolya=f minlen=0
	flaglongreads=f jni=f simd=f pigz=f unpigz=f qin=33 qout=33 ow=f -Xmx128m)
printf '%q ' "${normal[@]}" > "$out/normal.command"; printf '\n' >> "$out/normal.command"
printf '%q ' "${invalid[@]}" > "$out/invalid.command"; printf '\n' >> "$out/invalid.command"
normal_code=0; "${normal[@]}" > "$out/normal.stdout" 2> "$out/normal.stderr" || normal_code=$?
invalid_code=0; "${invalid[@]}" > "$out/invalid.stdout" 2> "$out/invalid.stderr" || invalid_code=$?
printf '%s\n' "$normal_code" > "$out/normal.exit"
printf '%s\n' "$invalid_code" > "$out/invalid.exit"
failures=0
printf 'case\texpected\tobserved\tresult\n'
if [[ $normal_code = 0 ]] && cmp -s "$out/normal.expected.fq" "$out/normal.fq"; then
	printf 'normal_output\t%dreads/%dbases/exactbytes\texact\tPASS\n' "$((wells*3))" "$((wells*12))"
else
	printf 'normal_output\texactbytes\texit%s_or_byte_difference\tFAIL\n' "$normal_code"
	failures=$((failures+1))
fi
if [[ $phase = baseline ]]; then
	printf 'invalid_worker_baseline\texit_and_diagnostics_recorded\texit%s\tCAPTURED\n' "$invalid_code"
elif [[ $invalid_code = 1 ]] && grep -Fq "$diagnostic" "$out/invalid.stderr" &&
	grep -Fq "$class worker 0 failed; aborting blocked input/output threads." "$out/invalid.stderr"; then
	printf 'invalid_worker\texit1/original_exception/fatal_worker_message\tmatched\tPASS\n'
else
	printf 'invalid_worker\texit1/original_exception/fatal_worker_message\texit%s_or_missing_diagnostic\tFAIL\n' "$invalid_code"
	failures=$((failures+1))
fi
[[ $failures = 0 ]]
