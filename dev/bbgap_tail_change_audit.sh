#!/usr/bin/env bash
# Joins changed paired-stop rows to normalized reference/baseline/candidate rows.
set -euo pipefail
runtime=${BBGAP_RUNTIME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
exec java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" prot.HbmTailChangeAudit "$@"
