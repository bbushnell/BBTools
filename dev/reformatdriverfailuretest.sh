#!/usr/bin/env bash
# Modern Reformat driver worker-failure control: exact normal output plus
# per-read tag parsing failure in a worker thread.
set -eo pipefail

SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)

out=
phase=fixed
bbtools=$DIR
driver=streamer
workers=2

for arg in "$@"; do
	case "$arg" in
	outdir=*) out=${arg#outdir=};;
	bbtools=*) bbtools=${arg#bbtools=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	driver=streamer|driver=reformat2) driver=${arg#driver=};;
	workers=1|workers=2) workers=${arg#workers=};;
	*) echo 'Use outdir=new-child-of-existing-directory [phase=baseline|fixed] [driver=streamer|reformat2] [workers=1|2] [bbtools=absolute-root]' >&2; exit 1;;
	esac
done

[[ $bbtools = /* && -f "$bbtools/reformat.sh" && -d "$bbtools/current" ]] || {
	echo 'bbtools must be an absolute root containing reformat.sh and current/' >&2
	exit 1
}
[[ -n $out && ! -e $out && -d ${out%/*} ]] || {
	echo 'outdir must be a new child of an existing directory' >&2
	exit 1
}

mkdir -- "$out"
out=$(cd "$out" && pwd)

normal_records=64
invalid_records=64
[[ $phase != baseline ]] || invalid_records=1

normal_in="$out/normal.fq"
invalid_in="$out/invalid.fq"
expected="$out/normal.expected.fq"
: > "$normal_in"
: > "$invalid_in"
: > "$expected"

for ((i=0; i<normal_records; i++)); do
	printf -v id 'read%03d' "$i"
	printf '@%s depth:2,\nACGT\n+\nIIII\n' "$id" >> "$normal_in"
	printf '@%s depth:2,\nACGT\n+\nIIII\n' "$id" >> "$expected"
done

for ((i=0; i<invalid_records; i++)); do
	printf -v id 'bad%03d' "$i"
	if [[ $i = 0 ]]; then
		printf '@%s depth:bad,\nACGT\n+\nIIII\n' "$id" >> "$invalid_in"
	else
		printf '@%s depth:2,\nACGT\n+\nIIII\n' "$id" >> "$invalid_in"
	fi
done

class=ReformatStreamer
if [[ $driver = reformat2 ]]; then
	class=Reformat2
	. "$bbtools/javasetup.sh"
	. "$bbtools/memdetect.sh"
	parseJavaArgs "--xmx=128m" "--xms=32m" "--mode=fixed" simd=f
	setEnvironment
	java_flags=($EA $EOOM $SIMD $XMX $XMS)
	base=(java "${java_flags[@]}" -cp "$bbtools/current/" jgi.Reformat2)
else
	base=(bash "$bbtools/reformat.sh" -Xmx128m -Xms32m)
fi

common=("wt=$workers" t=2 tin=1 tout=1 readbufferlength=1 readbuffers=2
	int=f qin=33 qout=33 minlength=0 simd=f
	tag=depth: delimiter=, minvalue=1
	pigz=f unpigz=f)
normal=("${base[@]}" "in=$normal_in" "out=$out/normal.out.fq" "${common[@]}" ow=t)
invalid=("${base[@]}" "in=$invalid_in" "out=$out/invalid.out.fq" "${common[@]}" ow=t)

printf '%q ' "${normal[@]}" > "$out/normal.command"
printf '\n' >> "$out/normal.command"
printf '%q ' "${invalid[@]}" > "$out/invalid.command"
printf '\n' >> "$out/invalid.command"

normal_code=0
"${normal[@]}" > "$out/normal.stdout" 2> "$out/normal.stderr" || normal_code=$?
invalid_code=0
"${invalid[@]}" > "$out/invalid.stdout" 2> "$out/invalid.stderr" || invalid_code=$?

printf '%s\n' "$normal_code" > "$out/normal.exit"
printf '%s\n' "$invalid_code" > "$out/invalid.exit"

failures=0
printf 'case\tdriver\texpected\tobserved\tresult\n'

if [[ $normal_code = 0 ]] && cmp -s "$expected" "$out/normal.out.fq"; then
	printf 'normal_output\t%s\t%dreads/%dbases/exactbytes\texact\tPASS\n' "$driver" "$normal_records" "$((normal_records*4))"
else
	printf 'normal_output\t%s\texactbytes\texit%s_or_byte_difference\tFAIL\n' "$driver" "$normal_code"
	failures=$((failures+1))
fi

fatal_message="$class worker [01] failed; aborting blocked input/output threads\\."
if grep -Fq 'AssertionError: 50 = b' "$out/invalid.stderr" &&
	grep -Fq 'at parse.Parse.parseDouble' "$out/invalid.stderr" &&
	grep -Fq 'at processor.ReformatProcessor.passesTagFilter' "$out/invalid.stderr"; then
	original_error=present
else
	original_error=missing
fi

if [[ $phase = baseline ]]; then
	if [[ $invalid_code = 1 && $original_error = present ]]; then
		printf 'invalid_worker_baseline\t%s\texit1_and_original_numeric_error\texit1/original_numeric_error=present\tPASS\n' "$driver"
	else
		printf 'invalid_worker_baseline\t%s\texit1_and_original_numeric_error\texit%s/original_numeric_error=%s\tFAIL\n' \
			"$driver" "$invalid_code" "$original_error"
		failures=$((failures+1))
	fi
elif [[ $invalid_code = 1 && $original_error = present ]] &&
	grep -Eq "$fatal_message" "$out/invalid.stderr"; then
	printf 'invalid_worker\t%s\texit1/original_numeric_error/fatal_worker_message\tmatched\tPASS\n' "$driver"
else
	printf 'invalid_worker\t%s\texit1/original_numeric_error/fatal_worker_message\texit%s/original_numeric_error=%s/fatal_message_missing\tFAIL\n' \
		"$driver" "$invalid_code" "$original_error"
	failures=$((failures+1))
fi

[[ $failures = 0 ]]
