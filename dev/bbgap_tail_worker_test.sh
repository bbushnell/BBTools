#!/usr/bin/env bash
# Real-library tail Worker parity checks. Run on Dori with a frozen runtime;
# ordinary no-argument search fixtures remain in bbgap_tail_profile_test.sh.
set -euo pipefail
runtime=${BBGAP_RUNTIME:?Set BBGAP_RUNTIME to the complete frozen BBTools runtime}
[[ -f $runtime/current/test/prot/ProteinTailSearchTest.class ]]
[[ $# == 1 && $1 == models=* && ${1#models=} != '' ]] || {
 echo 'Usage: BBGAP_RUNTIME=/path/to/runtime bash bbgap_tail_worker_test.sh models=/pack/models[,/other/models]' >&2
 exit 2
}
exec java -ea -Xmx8g -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" test.prot.ProteinTailSearchTest "$@"
