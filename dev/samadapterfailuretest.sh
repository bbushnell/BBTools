#!/usr/bin/env bash
# Run each normal/failure adapter case in a fresh JVM; caller checks failure exit and diagnostics.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ ${1:-} = suite ]]; then
	[[ $# = 2 && $2 = outdir=* ]] || { echo 'Use suite outdir=new-directory' >&2; exit 1; }
	test_dir=${2#outdir=}
	[[ -n $test_dir && ! -e $test_dir && -d ${test_dir%/*} ]] || exit 1
	mkdir -- "$test_dir"
	bash "$DIR/dev/samadapterfailuretest.sh" "out=$test_dir/normal.sam" failure=f \
		> "$test_dir/normal.stdout" 2> "$test_dir/normal.stderr"
	grep -Fxq 'SAM_ADAPTER_NORMAL_PASS records=1' "$test_dir/normal.stdout"
	failure_code=0
	bash "$DIR/dev/samadapterfailuretest.sh" "out=$test_dir/failure.sam" failure=t \
		> "$test_dir/failure.stdout" 2> "$test_dir/failure.stderr" || failure_code=$?
	printf '%s\n' "$failure_code" > "$test_dir/failure.exit"
	[[ $failure_code = 1 ]]
	grep -Fq 'Fatal SAM/BAM writer failure; output may be incomplete.' "$test_dir/failure.stderr"
	grep -Fq 'AssertionError: Sam format cannot be used as input' "$test_dir/failure.stderr"
	grep -Fq 'at stream.SamLine.<init>' "$test_dir/failure.stderr"
	grep -Fq 'at stream.ReadStreamSamWriter.processJobs' "$test_dir/failure.stderr"
	echo 'SAM_ADAPTER_LIFECYCLE_PASS normal_records=1 failure_exit=1'
	exit 0
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.stream.SamAdapterFailureTest "$@"
