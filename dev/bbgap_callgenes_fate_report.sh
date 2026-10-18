#!/usr/bin/env bash
# Summarizes CallGenes hbmdiag supported-stop fate rows; no caller execution.
set -euo pipefail
runtime=${BBGAP_RUNTIME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
exec java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$runtime/current" prot.ProteinSupportFateReport "$@"
