#!/usr/bin/env bash
set -euo pipefail
# families= defaults to RF00174,RF01482,RF01689. Other families require one
# matching loci= ledger per Rfam ID; sampling/masking never establishes absence.
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=2 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.CobalaminShredSampler "$@"
