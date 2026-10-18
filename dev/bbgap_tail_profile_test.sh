#!/usr/bin/env bash
# Literal profile-policy checks plus unchanged alignment/score regressions.
set -euo pipefail
build=${1:?frozen build}
[[ $build == /tmp/bbgap-callgenes-* && -f $build/PASS ]]
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bash "$repo/dev/bbgap_hbm_tail_test.sh" "$build"
for suite in test.prot.ProteinRowDotTest test.prot.HbmTailProfileTest test.prot.ProteinTailSearchTest test.prok.ProteinSupportTest test.prok.CdsTranslatorTest test.prok.ProteinGuidedModelTest; do
 java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
  -cp "$build/runtime/current" "$suite" > "$build/$suite.stdout" 2> "$build/$suite.stderr"
 cat "$build/$suite.stderr"
done
