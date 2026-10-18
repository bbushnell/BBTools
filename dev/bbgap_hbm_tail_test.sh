#!/usr/bin/env bash
# Runs exact alignment oracles in a frozen, complete-package build.
set -euo pipefail
build=${1:?frozen bbgap-callgenes build root}
[[ $build == /tmp/bbgap-callgenes-* && -f $build/PASS ]]
for test in prot.HbmPositionTest test.prot.HbmTailTest test.prot.HbmTailGateTest test.prot.HbmCompetitiveGateTest prot.FamilyShortlistScoresTest; do
	log=${test##*.}
	java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
		-cp "$build/runtime/current" "$test" > "$build/$log.stdout" 2> "$build/$log.stderr"
	cat "$build/$log.stderr"
done
