#!/usr/bin/env bash
# Joins rejected original tail rows to raw reference and selected GFF stops.
# Diagnostic only: raw reference-stop coordinates are not normalized truth.
set -euo pipefail
runtime=${BBGAP_RUNTIME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
exec java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" prot.HbmTailLossAudit "$@"
