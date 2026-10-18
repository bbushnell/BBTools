#!/usr/bin/env bash
# Early tail-cap CLI validation against an immutable full build; no HBM resources.
# Usage: bash dev/bbgap_callgenes_tailcap_test.sh /tmp/bbgap-callgenes-BUILD
# A deliberately invalid training selector is the immediately following guard:
# valid caps must reach it; invalid caps must fail first. Neither arm calls genes
# or loads default references. Actual cap geometry/GFF behavior needs Java tests
# and real-reference integration. Explicit trnalib/trnamodel loads are out of scope.
set -euo pipefail
build=${1:?frozen full CallGenes build required}
fail(){ echo "TAILCAP_CLI_FAIL: $*" >&2; exit 1; }
[[ $# == 1 && $build =~ ^/tmp/bbgap-callgenes-[a-zA-Z0-9_.-]+$ ]] || fail 'Expected one safe frozen-build path'
[[ -f $build/PASS && -f $build/source_revision.txt && -f $build/runtime/callgenes.sh \
 && -f $build/runtime/current/prok/CallGenes.class ]] || fail 'Incomplete frozen build'
grep -Fxq CALLGENES_SUPPORT_BUILD_PASS "$build/PASS" || fail 'Missing full-build PASS marker'
command -v timeout >/dev/null || fail 'GNU timeout required for the early-exit diagnostic watchdog'
root=$(mktemp -d "$build/tailcap_cli.XXXXXX")
printf 'Tail-cap CLI evidence: %s\n' "$root"
cp "$build/source_revision.txt" "$root/source_revision.txt"
cp "${BASH_SOURCE[0]}" "$root/recipe.sh"
printf '>fixture\nATGAAATAA\n' > "$root/input.fa"
cp "$root/input.fa" "$root/original.fa"
hostname > "$root/host"
date -u +%FT%TZ > "$root/started.utc"
command -v java > "$root/java.path"
java -version > "$root/java.version" 2>&1
missing=$root/missing_hbm
[[ ! -e $missing && ! -L $missing ]] || fail 'Missing-reference sentinel exists'
printf 'case\texit\tcheck\n' > "$root/results.tsv"
checks=0
trap 'rc=$?; printf "%s\n" "$rc" > "$root/exit"' EXIT

expect_failure(){
	local name=$1 kind=$2 annotation=$3 rule=$4 expected
	shift 4
	local dir=$root/$name rc=0 suffix
	case $kind in
		cap) expected='java[.]lang[.]IllegalArgumentException: Invalid tail length, terminal-support, identity or coverage thresholds' ;;
		parsed) expected='java[.]lang[.]IllegalArgumentException: hbmtrainingset must be auto, supported or selected: tailcap_probe' ;;
		overflow) expected='java[.]lang[.]AssertionError: Value 2147483648 is out of range for integers:' ;;
		*) fail "Unknown test kind: $kind" ;;
	esac
	mkdir "$dir"
	local -a cmd=(bash "$build/runtime/callgenes.sh" "in=$root/input.fa"
		"out=$dir/calls.gff" "outa=$dir/proteins.fa" "stats=$dir/stats.txt"
		"hist=$dir/hist.txt" "hbmdiag=$dir/fates.tsv" "hbmreference=$missing"
		"annotation=$annotation" "hbmtailrule=$rule" passes=2 hbmtrainingset=tailcap_probe
		cdsonly=t tax=f t=1 overwrite=f -Xmx256m -ea "$@")
	# This watchdog bounds only an expected immediate CLI rejection, not a science job.
	printf '%q ' timeout --kill-after=5s 30s "${cmd[@]}" > "$dir/argv.txt"
	printf '\n' >> "$dir/argv.txt"
	timeout --kill-after=5s 30s "${cmd[@]}" </dev/null > "$dir/stdout" 2> "$dir/stderr" || rc=$?
	printf '%s\n' "$rc" > "$dir/exit.txt"
	[[ $rc -ne 0 && $rc -lt 124 ]] || fail "$name: success, timeout, signal or launcher failure (exit=$rc)"
	# Match an exception line, never the launcher's echoed argv containing a sentinel.
	grep -E '(^|[[:space:]])java[.]lang[.](IllegalArgumentException|AssertionError):' "$dir/stderr" \
		> "$dir/exceptions.txt" || fail "$name: missing expected exception; see $dir/stderr"
	grep -Eq "$expected" "$dir/exceptions.txt" || fail "$name: wrong rejection reason; see $dir/exceptions.txt"
	if grep -Eq 'Missing BBGAP reference:|Could not load annotation reference:|Protein support: mode=' "$dir/stderr"; then
		fail "$name: reached HBM loading instead of the early validation guard"
	fi
	for suffix in calls.gff proteins.fa stats.txt hist.txt fates.tsv; do
		[[ ! -e $dir/$suffix && ! -L $dir/$suffix ]] || fail "$name: opened output before rejection: $suffix"
	done
	printf '%s\t%s\t%s\n' "$name" "$rc" "$kind" >> "$root/results.tsv"
	checks=$((checks+1))
}

for annotation in core full; do
	for rule in geometry profile; do
		for cap in 30 200 300 400 500 600; do
			expect_failure "${annotation}_${rule}_${cap}" parsed "$annotation" "$rule" "hbmtailcap=$cap"
		done
	done
done
for annotation in off core full; do
	for cap in -1 0 29; do
		expect_failure "${annotation}_invalid_${cap}" cap "$annotation" geometry "hbmtailcap=$cap"
	done
done
expect_failure profile_invalid_29 cap core profile hbmtailcap=29
expect_failure empty_cap cap core geometry hbmtailcap=
expect_failure null_cap cap core geometry hbmtailcap=null
expect_failure integer_overflow overflow core geometry hbmtailcap=2147483648
expect_failure default_geometry parsed core geometry
expect_failure default_profile parsed core profile
expect_failure kmg_cap parsed core geometry hbmtailcap=0.2k
[[ $checks == 40 ]] || fail "Incomplete fixture matrix: $checks checks"
cmp "$root/input.fa" "$root/original.fa"
[[ ! -e $missing && ! -L $missing ]] || fail 'Missing-reference sentinel changed'
date -u +%FT%TZ > "$root/finished.utc"
printf 'CALLGENES_TAILCAP_CLI_PASS checks=%s invalid_caps=true allowed_caps_parse=true outputs_unopened=true scope=early_validation_only\n' "$checks" > "$root/PASS"
cat "$root/PASS"
