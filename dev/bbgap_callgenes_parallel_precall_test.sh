#!/usr/bin/env bash
# Bounded exact-parity checks; root prepares first-four-contig B4/A4 inputs with
# the existing reformat.sh and runs this on Dori against already deployed builds.
# No compilation, submission, data discovery, or modification of either runtime.
# Usage: selftest | failure DIAGNOSTIC_BUILD NEW_ROOT
#        run CANDIDATE_BUILD D86_BUILD CALLER_CONFIG A.fa B.fa NEW_ROOT
set -euo pipefail
export LC_ALL=C
fail(){ echo "PARALLEL_PRECALL_FAIL: $*" >&2; exit 1; }
record(){ printf '%q ' "$@"; printf '\n'; }

# Keep every pre-call field except input path and wall time, in emitted order.
precall(){
	local file=$1 line token sep count=0
	local -a fields
	while IFS= read -r line; do
		[[ $line == 'Protein-guided training: '* ]] || continue
		read -ra fields <<< "${line#Protein-guided training: }"
		sep=''
		for token in "${fields[@]}"; do
			case $token in input=*|seconds=*) continue ;; esac
			printf '%s%s' "$sep" "$token"; sep=' '
		done
		printf '\n'; count=$((count+1))
	done < "$file"
	[[ $count -gt 0 ]]
}

fates(){
	local dir=$1
	head -n 1 "$dir/fate.tsv" > "$dir/fate.header"
	grep -q $'^schema\tstop_id\t' "$dir/fate.header" || fail "Bad fate header: $dir"
	tail -n +2 "$dir/fate.tsv" | sort -t $'\t' -k2,2 > "$dir/fate.sorted"
	cut -f2 "$dir/fate.sorted" | uniq -d > "$dir/fate.duplicate_keys"
	[[ ! -s $dir/fate.duplicate_keys ]] || fail "Duplicate fate keys: $dir"
}

case ${1:-} in
selftest)
	[[ $# == 1 ]] || fail 'selftest accepts no other arguments'
	fixture=$(mktemp -d /tmp/bbgap-parallel-precall-selftest.XXXXXX)
	printf 'noise\nProtein-guided training: input=A contigs=4 bases=99 selected_supported=2 seconds=0.1 adapted=true training_set=supported\nProtein-guided training: input=B contigs=3 bases=88 seconds=9 adapted=false training_set=selected\n' > "$fixture/stderr"
	printf 'contigs=4 bases=99 selected_supported=2 adapted=true training_set=supported\ncontigs=3 bases=88 adapted=false training_set=selected\n' > "$fixture/expected"
	precall "$fixture/stderr" > "$fixture/actual"; cmp "$fixture/expected" "$fixture/actual"
	printf 'no training\n' > "$fixture/absent"
	if precall "$fixture/absent"; then fail 'Missing training block accepted'; fi
	printf 'schema\tstop_id\tfate\n1\tb\tselected\n1\ta\tdp_not_selected\n' > "$fixture/fate.tsv"
	fates "$fixture"
	printf '1\ta\tdp_not_selected\n1\tb\tselected\n' > "$fixture/fate.expected"
	cmp "$fixture/fate.expected" "$fixture/fate.sorted"
	printf '1\ta\tselected\n' >> "$fixture/fate.tsv"
	if (fates "$fixture") > "$fixture/duplicate.stdout" 2> "$fixture/duplicate.stderr"; then fail 'Duplicate fate key accepted'; fi
	printf 'PARALLEL_PRECALL_SELFTEST_PASS fixture=%s\n' "$fixture"
	exit 0 ;;
failure)
	[[ $# == 3 ]] || fail 'failure requires diagnostic build and fresh evidence root'
	build=$2 root=$3
	for path in "$build" "$root"; do
		[[ $path =~ ^/[A-Za-z0-9_./-]+$ && $path != */../* && $path != */./* ]] || fail "Unsafe path: $path"
	done
	grep -Fxq CALLGENES_SUPPORT_BUILD_PASS "$build/PASS"
	[[ -f $build/runtime/current/test/prok/ProteinTrainingFailureTest.class && ! -e $root && ! -L $root ]]
	mkdir "$root"
	printf '>worker_failure_fixture\nATGAAATAA\n' > "$root/input.fa"
	cp "$build/source_revision.txt" "$root/source_revision.txt"
	command=(nice -n10 java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector
		-cp "$build/runtime/current" test.prok.ProteinTrainingFailureTest "$root/input.fa")
	record timeout --kill-after=5s 30s "${command[@]}" > "$root/argv.txt"
	hostname > "$root/host"; date -u +%FT%TZ > "$root/started.utc"
	rc=0
	timeout --kill-after=5s 30s "${command[@]}" > "$root/stdout" 2> "$root/stderr" || rc=$?
	printf '%s\n' "$rc" > "$root/exit"
	[[ $rc == 1 ]] || fail "Expected whole-process halt1, not success or timeout; observed $rc"
	grep -Fxq 'PROTEIN_TRAINING_FAILURE_ARMED consumer_daemon=false' "$root/stderr"
	grep -Fq 'java.lang.IllegalStateException: PROTEIN_TRAINING_WORKER_SENTINEL' "$root/stderr"
	grep -Fq 'prok.CallGenes$ProteinTrainingThread.run(' "$root/stderr"
	grep -Fq 'Protein-guided input/calling failed for ' "$root/stderr"
	if grep -Eq '^Protein-guided training: |Training worker returned without halting' "$root/stderr"; then
		fail 'Failure path returned or reported completed training'
	fi
	date -u +%FT%TZ > "$root/finished.utc"
	printf 'PARALLEL_PRECALL_FAILURE_PASS sentinel=true actual_worker=true live_consumer=true halt_exit=1\n' > "$root/PASS"
	cat "$root/PASS"; exit 0 ;;
run) shift ;;
*) fail 'Expected selftest, failure or run' ;;
esac
[[ $# == 6 ]] || fail 'Expected candidate, baseline, config, A, B, fresh root'
candidate=$1 baseline=$2 config=$3 a=$4 b=$5 root=$6
for path in "$@"; do
	[[ $path =~ ^/[A-Za-z0-9_./-]+$ && $path != */../* && $path != */./* ]] || fail "Unsafe path: $path"
done
[[ ! -e $root && ! -L $root && -f $config && -f $a && -f $b && $a != "$b" ]] || fail 'Missing inputs or existing output root'
for build in "$candidate" "$baseline"; do
	grep -Fxq CALLGENES_SUPPORT_BUILD_PASS "$build/PASS"
	[[ -s $build/source_revision.txt && -f $build/runtime/current/prok/CallGenes.class && -f $build/runtime/callgenes.sh ]]
done
grep -Eq '^d86dd3f1[0-9a-f]{32}$' "$baseline/source_revision.txt" || fail 'Baseline must be the frozen training-set runtime'
[[ $(grep -c '^>' "$a") == 4 && $(grep -c '^>' "$b") == 4 ]] || fail 'A and B must each contain exactly four prepared FASTA records'
command -v timeout >/dev/null
mkdir -p "$root/a2"
cp "$a" "$root/a2/${a##*/}"; a2=$root/a2/${a##*/}; cmp "$a" "$a2"
cp "$config" "$root/caller.config"
cp "$candidate/source_revision.txt" "$root/candidate_revision.txt"
cp "$baseline/source_revision.txt" "$root/baseline_revision.txt"
record bash "${BASH_SOURCE[0]}" run "$@" > "$root/argv.txt"
hostname > "$root/host"; date -u +%FT%TZ > "$root/started.utc"
java -version > "$root/java.version" 2>&1
printf 'arm\texit\n' > "$root/exits.tsv"
trap 'rc=$?; printf "%s\n" "$rc" > "$root/exit"' EXIT

run_arm(){
	local name=$1 build=$2 labels=$3 input=$4 threads=$5 buffer=$6 limit=$7 diagnostic=$8 rc=0
	local dir=$root/$name
	mkdir "$dir"
	local -a command=(nice -n10 bash "$build/runtime/callgenes.sh" "config=$root/caller.config"
		"in=$input" "out=$dir/calls.gff" annotation=core passes=2 "hbmtrainingset=$labels"
		tax=f ordered=t overwrite=f "t=$threads" "maxreads=$limit" "readbufferlen=$buffer"
		hbmtailrule=geometry hbmboostmode=flat hbmprotect=auto hbmboost=25 hbmstopboost=25
		hbmsearch=fast hbmminseeds=1 hbmminseeddensity=0 hbmlookahead=4
		-Xmx12g -Xms32m -eoom -ea)
	if [[ $diagnostic == yes ]]; then command+=("hbmdiag=$dir/fate.tsv"); else command+=(hbmdiag=null); fi
	record "${command[@]}" > "$dir/argv.txt"
	/usr/bin/time -f 'elapsed_seconds\t%e\nuser_seconds\t%U\nsystem_seconds\t%S' -o "$dir/time.tsv" \
		"${command[@]}" > "$dir/stdout" 2> "$dir/stderr" || rc=$?
	printf '%s\n' "$rc" > "$dir/exit"; printf '%s\t%s\n' "$name" "$rc" >> "$root/exits.tsv"
	[[ $rc == 0 ]] || fail "$name exited $rc"
	[[ -s $dir/calls.gff ]] || fail "Missing GFF: $name"
	precall "$dir/stderr" > "$dir/precall.tsv"
	grep -E '^Protein (stop support|tail search work):' "$dir/stderr" > "$dir/final_counters.txt"
	[[ $(wc -l < "$dir/final_counters.txt") == 2 ]] || fail "Missing final counters: $name"
	if [[ $diagnostic == yes ]]; then fates "$dir"; fi
}

compare(){
	local left=$root/$1 right=$root/$2 diagnostic=$3 file
	for file in calls.gff precall.tsv final_counters.txt; do cmp "$left/$file" "$right/$file"; done
	if [[ $diagnostic == yes ]]; then
		cmp "$left/fate.header" "$right/fate.header"
		cmp "$left/fate.sorted" "$right/fate.sorted"
	fi
}

# Shared.java READ_BUFFER_LENGTH=200 and bufferLen() returns it. Pin that normal
# batch size explicitly, rather than inheriting a possibly overridden config value.
for labels in supported selected; do
	run_arm "${labels}_old_a" "$baseline" "$labels" "$a" 1 200 4 yes
	grep -Eq '(^| )adapted=true( |$)' "$root/${labels}_old_a/precall.tsv" || fail 'A must exercise actual model adaptation'
	grep -Eq '(^| )selected_supported=[1-9][0-9]*( |$)' "$root/${labels}_old_a/precall.tsv" || fail 'A must include supported final stops'
	run_arm "${labels}_new_a1" "$candidate" "$labels" "$a" 1 200 4 yes
	run_arm "${labels}_new_a4" "$candidate" "$labels" "$a" 4 200 4 yes
	run_arm "${labels}_new_batches" "$candidate" "$labels" "$a" 4 1 4 yes
	for name in new_a1 new_a4 new_batches; do compare "${labels}_old_a" "${labels}_$name" yes; done
	grep -Eq '^contigs=4 ' "$root/${labels}_new_batches/precall.tsv"
	run_arm "${labels}_old_limit" "$baseline" "$labels" "$a" 4 1 3 yes
	run_arm "${labels}_new_limit" "$candidate" "$labels" "$a" 4 1 3 yes
	compare "${labels}_old_limit" "${labels}_new_limit" yes
	grep -Eq '^contigs=3 ' "$root/${labels}_new_limit/precall.tsv"
	run_arm "${labels}_old_b" "$baseline" "$labels" "$b" 1 200 4 yes
	run_arm "${labels}_new_b" "$candidate" "$labels" "$b" 4 1 4 yes
	compare "${labels}_old_b" "${labels}_new_b" yes
	# Repeated A has the same stop keys; the intentional fate duplicate guard
	# prohibits hbmdiag here. Individual A/B runs above compare complete fate rows.
	run_arm "${labels}_old_aba" "$baseline" "$labels" "$a,$b,$a2" 4 1 4 no
	run_arm "${labels}_new_aba" "$candidate" "$labels" "$a,$b,$a2" 4 1 4 no
	compare "${labels}_old_aba" "${labels}_new_aba" no
	[[ $(wc -l < "$root/${labels}_new_aba/precall.tsv") == 3 ]]
	cmp <(sed -n '1p' "$root/${labels}_new_aba/precall.tsv") <(sed -n '3p' "$root/${labels}_new_aba/precall.tsv")
	cat "$root/${labels}_new_batches/precall.tsv" "$root/${labels}_new_b/precall.tsv" "$root/${labels}_new_batches/precall.tsv" > "$root/${labels}_standalone.precall.tsv"
	cmp "$root/${labels}_standalone.precall.tsv" "$root/${labels}_new_aba/precall.tsv"
	{ cat "$root/${labels}_new_batches/calls.gff"; sed '/^##/d' "$root/${labels}_new_b/calls.gff" "$root/${labels}_new_batches/calls.gff"; } > "$root/${labels}_standalone.gff"
	cmp "$root/${labels}_standalone.gff" "$root/${labels}_new_aba/calls.gff"
done
date -u +%FT%TZ > "$root/finished.utc"
printf 'PARALLEL_PRECALL_PARITY_PASS label_modes=2 success_arms=20 full_gff=true sorted_fates=true precall_counters=true final_counters=true maxreads=true aba=true worker_catch=UNTESTED\n' > "$root/PASS"
cat "$root/PASS"
