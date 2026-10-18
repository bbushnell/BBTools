#!/usr/bin/env bash
# Sayu: early annotated passes=2 CLI guards against an immutable full build.
# Usage: bash dev/bbgap_callgenes_guided_test.sh /tmp/bbgap-callgenes-BUILD
# No real HBM package, compilation, network access, or scheduler work. Retain
# fixtures/argv/exits/logs under the build; PASS covers guards, not adaptation.
set -euo pipefail
build=${1:?frozen full CallGenes build required}
fail(){ echo "GUIDED_CLI_FAIL: $*" >&2; exit 1; }
[[ $# == 1 ]] || fail 'Expected exactly one frozen-build argument'
# The existing launcher evaluates its command string; keep generated paths literal.
[[ $build =~ ^/tmp/bbgap-callgenes-[a-zA-Z0-9_.-]+$ ]] || fail "Unsafe build path: $build"
[[ -f $build/PASS && -f $build/source_revision.txt && -f $build/runtime/callgenes.sh \
 && -f $build/runtime/current/prok/CallGenes.class ]] || fail 'Incomplete frozen build'
grep -Fxq CALLGENES_SUPPORT_BUILD_PASS "$build/PASS" || fail 'Missing full-build PASS marker'
command -v timeout >/dev/null || fail 'GNU timeout is required to detect blocked input guards'
root=$(mktemp -d "$build/guided_cli.XXXXXX")
printf 'Guided CLI evidence: %s\n' "$root"
cp "$build/source_revision.txt" "$root/source_revision.txt"
printf '>fixture\nATGAAATAA\n' > "$root/input.fa"
printf '@fixture\nATGAAATAA\n+\nIIIIIIIII\n' > "$root/input.fq"
printf '##gff-version 3\nfixture\ttest\tCDS\t1\t9\t.\t+\t0\tID=fixture\n' > "$root/input.gff"
cp "$root/input.fa" "$root/original.fa"
cp "$root/input.fq" "$root/original.fq"
cp "$root/input.gff" "$root/original.gff"
gzip -c "$root/input.fa" > "$root/input.fa.gz"
ln -s "$root/input.fa" "$root/linked.fa"
mkfifo "$root/pipe.fa"
missing=$root/missing_hbm
[[ ! -e $missing ]] || fail 'Missing-reference sentinel unexpectedly exists'
printf 'case\texit\tcheck\n' > "$root/results.tsv"
checks=0

# Match only an exception message, never the launcher's echoed argv. The absent
# HBM root distinguishes an early policy rejection from model-loading failure.
expect_failure(){
	local name=$1 expected=$2 kind=$3 input=$4
	shift 4
	local dir=$root/$name rc=0 suffix
	mkdir "$dir"
	local -a cmd=(bash "$build/runtime/callgenes.sh" "in=$input"
		"out=$dir/calls.gff" "outa=$dir/proteins.fa" "stats=$dir/stats.txt"
		"hist=$dir/hist.txt" "hbmdiag=$dir/fates.tsv" "hbmreference=$missing"
		cdsonly=t tax=f t=1 overwrite=f -Xmx256m -ea "$@")
	printf '%q ' timeout --kill-after=5s 30s "${cmd[@]}" > "$dir/argv.txt"
	printf '\n' >> "$dir/argv.txt"
	timeout --kill-after=5s 30s "${cmd[@]}" </dev/null > "$dir/stdout" 2> "$dir/stderr" || rc=$?
	printf '%s\n' "$rc" > "$dir/exit.txt"
	[[ $rc -ne 0 && $rc -lt 124 ]] || fail "$name: success, timeout, signal or launcher failure (exit=$rc); see $dir"
	grep -E '(^|[[:space:]])java\.lang\.IllegalArgumentException:' "$dir/stderr" \
		> "$dir/exceptions.txt" || fail "$name: no argument exception; see $dir/stderr"
	grep -Eq "$expected" "$dir/exceptions.txt" || fail "$name: wrong rejection reason; see $dir/exceptions.txt"
	if [[ $kind == guard ]]; then
		if grep -Eq 'Missing BBGAP reference:|Could not load annotation reference:|Protein support: mode=' "$dir/stderr"; then
			fail "$name: reached HBM loading instead of rejecting the input policy"
		fi
	elif [[ $kind != loader ]]; then
		fail "$name: unknown test kind $kind"
	fi
	for suffix in calls.gff proteins.fa stats.txt hist.txt fates.tsv; do
		[[ ! -e $dir/$suffix && ! -L $dir/$suffix ]] || fail "$name: output opened before rejection: $suffix"
	done
	printf '%s\t%s\t%s\n' "$name" "$rc" "$kind" >> "$root/results.tsv"
	checks=$((checks+1))
}

for mode in core full; do
	# Legal combinations must get PAST the new guards, reaching only the missing
	# reference sentinel. This is not successful biological calling validation.
	for pass in 1 2; do
		expect_failure "${mode}_allowed_${pass}" 'Missing BBGAP reference:' loader "$root/input.fa" "annotation=$mode" "passes=$pass"
	done
	expect_failure "${mode}_gzip" 'Missing BBGAP reference:' loader "$root/input.fa.gz" "annotation=$mode" passes=2
	expect_failure "${mode}_symlink" 'Missing BBGAP reference:' loader "$root/linked.fa" "annotation=$mode" passes=2
	for pass in -1 0 3; do
		expect_failure "${mode}_passes_${pass}" 'passes=1.*2' guard "$root/input.fa" "annotation=$mode" "passes=$pass"
	done
	for pass in 1 2; do
		for operation in recode retranslate; do
			expect_failure "${mode}_${operation}_${pass}" 'translation mode' guard "$root/input.fa" "annotation=$mode" "passes=$pass" "$operation=t"
		done
	done
	expect_failure "${mode}_ingff" 'ingff' guard "$root/input.fa" "annotation=$mode" passes=2 "ingff=$root/input.gff"
	for option in interleaved merge ecco; do
		expect_failure "${mode}_${option}" 'interleav|merge|ecco' guard "$root/input.fa" "annotation=$mode" passes=2 "$option=t"
	done
	expect_failure "${mode}_stdin" 'regular|replayable' guard stdin.fa "annotation=$mode" passes=2
	expect_failure "${mode}_fifo" 'regular|replayable' guard "$root/pipe.fa" "annotation=$mode" passes=2
	expect_failure "${mode}_fastq" 'FASTA' guard "$root/input.fq" "annotation=$mode" passes=2
	# Check every input, not merely the first member of a multi-genome request.
	expect_failure "${mode}_second_fastq" 'FASTA' guard "$root/input.fa,$root/input.fq" "annotation=$mode" passes=2
	# Training-set policy must reject unsupported values/combinations before HBM loading.
	expect_failure "${mode}_training_invalid" 'hbmtrainingset.*supported.*selected' guard "$root/input.fa" "annotation=$mode" passes=2 hbmtrainingset=invalid
	expect_failure "${mode}_training_selected_p1" 'hbmtrainingset=selected.*passes=2' guard "$root/input.fa" "annotation=$mode" passes=1 hbmtrainingset=selected
	expect_failure "${mode}_training_selected_p2" 'Missing BBGAP reference:' loader "$root/input.fa" "annotation=$mode" passes=2 hbmtrainingset=selected
done
# Disable the otherwise-invalid off-mode fate log, so it cannot mask this guard.
expect_failure off_training_selected 'hbmtrainingset=selected.*annotation' guard "$root/input.fa" annotation=off passes=2 hbmtrainingset=selected hbmdiag=null
cmp "$root/input.fa" "$root/original.fa"
cmp "$root/input.fq" "$root/original.fq"
cmp "$root/input.gff" "$root/original.gff"
printf 'CALLGENES_GUIDED_CLI_PASS checks=%s forbidden_inputs=true unsupported_passes=true outputs_unopened=true\n' "$checks" > "$root/PASS"
cat "$root/PASS"
