#!/usr/bin/env bash
# Bounded GenerateMultiChrom arity and MergeReadsAndGenome writer-error fixture.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
outdir=
phase=baseline
for arg in "$@"; do
	case "$arg" in
	bbtools=*) DIR=${arg#bbtools=};;
	outdir=*) outdir=${arg#outdir=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	*) echo 'Use outdir=fresh-absolute-directory [phase=baseline|fixed] [bbtools=absolute-root]' >&2; exit 1;;
	esac
done
[[ $DIR = /* && -f $DIR/javasetup.sh && -d $DIR/current && -n $outdir && $outdir = /* && ! -e $outdir && -d ${outdir%/*} ]]
mkdir -- "$outdir"
outdir=$(cd "$outdir" && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment

bases=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
printf '>r0\n%s\n' "$bases" > "$outdir/input.fa"
printf '>0\n%s\n' "$bases" > "$outdir/expected.fa"
failures=0

run_case(){
	local name=$1
	local expected_exit=$2
	local required_stderr=$3
	shift 3
	local caseout=$outdir/$name
	mkdir -- "$caseout"
	printf '%q ' "$@" > "$caseout/command.txt"
	printf '\n' >> "$caseout/command.txt"
	local exitcode=0
	"$@" > "$caseout/stdout.txt" 2> "$caseout/stderr.txt" || exitcode=$?
	printf '%s\n' "$exitcode" > "$caseout/exit.txt"
	local result=PASS
	[[ $exitcode == "$expected_exit" ]] || result=FAIL
	[[ -z $required_stderr ]] || grep -q -- "$required_stderr" "$caseout/stderr.txt" || result=FAIL
	printf '%s\t%s\t%s\t%s\n' "$name" "$expected_exit" "$exitcode" "$result"
	[[ $result = PASS ]] || failures=$((failures+1))
}

expected_arity=ArrayIndexOutOfBoundsException
expected_full_exit=0
if [[ $phase = fixed ]]; then
	expected_arity=IllegalArgumentException
	expected_full_exit=1
fi

printf 'case\texpected_exit\tobserved_exit\tresult\n'
run_case generatemultichrom-arity 1 "$expected_arity" \
	java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
	pacbio.GenerateMultiChrom
run_case merge-normal 0 '' \
	java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
	pacbio.MergeReadsAndGenome in="$outdir/input.fa" out="$outdir/normal.fa" ow=t
if ! cmp -s "$outdir/expected.fa" "$outdir/normal.fa"; then
	printf 'merge-normal-output\tcmp\tmismatch\tFAIL\n'
	failures=$((failures+1))
fi
run_case merge-dev-full "$expected_full_exit" 'No space left on device' \
	java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
	pacbio.MergeReadsAndGenome in="$outdir/input.fa" out=/dev/full ow=t

exit "$failures"
