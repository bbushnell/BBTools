#!/usr/bin/env bash
# Bounded NN dependency checks for an immutable, full CallGenes build.
# Usage: bash dev/bbgap_callgenes_dependency_test.sh /tmp/bbgap-callgenes-BUILD
# Compiles ALL frozen current/test/ml sources, but runs ONLY the existing default
# normalization and weight-bits analytic fixtures through their public launchers.
# Evidence/classes live in a fresh /tmp directory; the supplied build is read-only.
# No biological calling, new oracle, mixed-mode experiment or CDS-NN concurrency
# claim. A future public handoff can place this helper in BBTools/dev/.
set -euo pipefail
fail(){ echo "CALLGENES_DEPENDENCY_FAIL: $*" >&2; exit 1; }
[[ $# == 1 ]] || fail 'Expected one immutable full-build root'
build=$1
[[ $build =~ ^/tmp/bbgap-callgenes-[A-Za-z0-9_.-]+$ && -d $build && ! -L $build ]] \
	|| fail 'Expected a literal /tmp/bbgap-callgenes-BUILD directory'
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
runtime=$build/runtime
grep -Fxq CALLGENES_SUPPORT_BUILD_PASS "$build/PASS"
revision=$(<"$build/source_revision.txt")
[[ $revision =~ ^[0-9a-f]{40}$ ]] || fail 'Missing full source revision'
[[ $(git -C "$repo" rev-parse "$revision^{commit}") == "$revision" ]] || fail 'Unknown build revision'
for file in current/prok/CallGenes.class current/ml/CellNet.class current/ml/Cell.class \
	current/test/ml/CellNetInputNormalizationTest.java current/test/ml/CellNetWeightBitsTest.java \
	javasetup.sh memdetect.sh; do
	[[ -f $runtime/$file ]] || fail "Incomplete frozen runtime: $file"
done
shopt -s nullglob
sources=("$runtime"/current/test/ml/*.java)
[[ ${#sources[@]} -gt 0 ]] || fail 'Frozen ML test package is empty'
for tool in java javac jar git bash; do command -v "$tool" >/dev/null || fail "Missing $tool"; done
java_bin=$(readlink -f "$(command -v java)")
javac_bin=$(readlink -f "$(command -v javac)")
jar_bin=$(readlink -f "$(command -v jar)")
[[ ${java_bin%/*} == "${javac_bin%/*}" && ${java_bin%/*} == "${jar_bin%/*}" ]] \
	|| fail 'java, javac and jar must belong to one JDK'

root=$(mktemp -d /tmp/bbgap-callgenes-dependency.XXXXXX)
printf 'Dependency evidence: %s\n' "$root"
trap 'rc=$?; printf "%s\n" "$rc" > "$root/exit"' EXIT
mkdir -p "$root/classes" "$root/runtime/dev"
cp "$build/source_revision.txt" "$root/source_revision.txt"
cp "${BASH_SOURCE[0]}" "$root/recipe.sh"
printf '%s\n' "$build" > "$root/build_root.txt"
printf '%q ' bash "${BASH_SOURCE[0]}" "$build" > "$root/argv.txt"
printf '\n' >> "$root/argv.txt"
printf '%s\n' "${sources[@]}" > "$root/test_sources.txt"
printf 'java\t%s\njavac\t%s\njar\t%s\n' "$java_bin" "$javac_bin" "$jar_bin" > "$root/tool_paths.tsv"
hostname > "$root/host"
date -u +%FT%TZ > "$root/started.utc"
printf 'step\texit\n' > "$root/exits.tsv"

# Both analytic main methods require zero application arguments. In particular,
# weightbits forwards t=1 to its fixture, which rejects it. ActiveProcessorCount=1
# supplies the one-CPU default without changing the public launchers or tests.
unset JDK_JAVA_OPTIONS _JAVA_OPTIONS CLASSPATH
export JAVA_TOOL_OPTIONS='-XX:ActiveProcessorCount=1'
printf '%s\n' 'JAVA_TOOL_OPTIONS=-XX:ActiveProcessorCount=1' \
	'JDK_JAVA_OPTIONS/_JAVA_OPTIONS/CLASSPATH unset; analytic mains receive zero args' \
	> "$root/jvm_environment.txt"
run(){
	local name=$1 rc=0; shift
	mkdir "$root/$name"
	printf '%q ' "$@" > "$root/$name/argv.txt"
	printf '\n' >> "$root/$name/argv.txt"
	date -u +%FT%TZ > "$root/$name/started.utc"
	"$@" > "$root/$name/stdout" 2> "$root/$name/stderr" || rc=$?
	printf '%s\n' "$rc" > "$root/$name/exit"
	date -u +%FT%TZ > "$root/$name/finished.utc"
	printf '%s\t%s\n' "$name" "$rc" >> "$root/exits.tsv"
	[[ $rc == 0 ]] || fail "$name exited $rc; see $root/$name/stderr"
}
run java_version "$java_bin" -version
run javac_version "$javac_bin" -version
run jar_version "$jar_bin" --version

# The ordinary full-build archive omits dev/. Recover only these two public
# launchers from its exact revision, not the evolving working tree.
for name in cellnetinputnormalizationtest cellnetweightbits; do
	run "source_$name" git -C "$repo" show "$revision:dev/$name.sh"
	cp "$root/source_$name/stdout" "$root/runtime/dev/$name.sh"
done
cp "$runtime/javasetup.sh" "$runtime/memdetect.sh" "$root/runtime/"

# Compile the entire test package without exclusions. Any needed source dependency
# comes from the SAME frozen revision, and every emitted class stays diagnostic.
run compile "$javac_bin" -J-Xmx1g -J-XX:ActiveProcessorCount=1 \
	--add-modules jdk.incubator.vector -cp "$runtime/current" \
	-sourcepath "$runtime/current" -d "$root/classes" "${sources[@]}"
[[ -f $root/classes/test/ml/CellNetInputNormalizationTest.class \
	&& -f $root/classes/test/ml/CellNetWeightBitsTest.class ]] || fail 'Missing diagnostic test classes'

# Public launchers support a jar classpath. A manifest-only jar adds diagnostic
# classes ahead of frozen dependencies without copying or overwriting the runtime.
# jar writes continuation lines for the long Class-Path header.
printf 'Manifest-Version: 1.0\nClass-Path: file:%s/classes/ file:%s/current/\n\n' \
	"$root" "$runtime" > "$root/classpath.mf"
run classpath "$jar_bin" cfm "$root/runtime/bbtools.jar" "$root/classpath.mf"
run normalization bash "$root/runtime/dev/cellnetinputnormalizationtest.sh" -Xmx256m -Xms32m -ea
grep -Fxq 'CellNetInputNormalizationTest PASS' "$root/normalization/stdout" \
	|| fail 'Normalization exited zero without its default-fixture PASS'
run weightbits bash "$root/runtime/dev/cellnetweightbits.sh" selftest=t -Xmx256m -Xms32m -ea
grep -q '^CellNetWeightBitsTest PASS ' "$root/weightbits/stdout" \
	|| fail 'Weight-bits exited zero without its fixture PASS'
date -u +%FT%TZ > "$root/finished.utc"
printf '%s\n' 'CALLGENES_DEPENDENCY_TEST_PASS full_test_ml_compile=true normalization=true weightbits=true available_cpus=1 biological_calling=UNTESTED optional_cds_nn_concurrency=UNTESTED' > "$root/PASS"
cat "$root/PASS"
