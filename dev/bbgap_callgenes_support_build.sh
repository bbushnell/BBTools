#!/usr/bin/env bash
# Builds an immutable CallGenes/AnnotateGenes runtime from a committed revision.
# Complete prok/prot packages, native JVM and Java 8 API compatibility.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
revision=${1:?source revision}
build=${2:?fresh absolute build directory}
[[ $build == /tmp/bbgap-callgenes-* && ! -e $build ]]
mkdir -p "$build/runtime" "$build/classes" "$build/java8"
git -C "$repo" rev-parse "$revision^{commit}" > "$build/source_revision.txt"
git -C "$repo" archive "$revision" current resources networks callgenes.sh annotategenes.sh javasetup.sh memdetect.sh \
	| tar --exclude='*.class' -xf - -C "$build/runtime"
javac -J-Xmx2g -J-XX:ActiveProcessorCount=2 --add-modules jdk.incubator.vector \
	-cp "$build/classes" -sourcepath "$build/runtime/current" -d "$build/classes" \
	"$build"/runtime/current/prot/*.java "$build"/runtime/current/prok/*.java \
	"$build"/runtime/current/test/prot/*.java "$build"/runtime/current/test/prok/*.java \
	> "$build/compile.stdout" 2> "$build/compile.stderr"
javac -J-Xmx2g -J-XX:ActiveProcessorCount=2 --release 8 \
	-cp "$build/classes" -sourcepath '' -d "$build/java8" \
	"$build"/runtime/current/prot/*.java "$build"/runtime/current/prok/*.java \
	"$build"/runtime/current/test/prot/*.java "$build"/runtime/current/test/prok/*.java \
	> "$build/java8.stdout" 2> "$build/java8.stderr"
cp -a "$build/classes/." "$build/runtime/current/"
printf 'CALLGENES_SUPPORT_BUILD_PASS\n' > "$build/PASS"
cat "$build/compile.stderr" "$build/java8.stderr" "$build/PASS"
