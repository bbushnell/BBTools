#!/usr/bin/env bash
# Summarizes stored tail diagnostic rows; no alignment or caller-policy change.
set -euo pipefail
runtime=${BBGAP_RUNTIME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
exec java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" prot.HbmTailPolicySummary "$@"
