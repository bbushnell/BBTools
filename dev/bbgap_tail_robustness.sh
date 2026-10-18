#!/usr/bin/env bash
# Direct known-family suffix perturbations; synthetic robustness, not accuracy.
# BBGAP_RUNTIME may point to an immutable complete build; normal CLI is flag=value.
set -euo pipefail
runtime=${BBGAP_RUNTIME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
exec java -ea -Xmx4g -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" prot.HbmTailRobustness "$@"
