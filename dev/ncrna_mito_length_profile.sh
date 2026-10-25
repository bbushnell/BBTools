#!/usr/bin/env bash
# Derive full-covered-member length bounds from an exact saved matrix; no alignment.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=4 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" prok.NcrnaMitoLengthProfile "$@"
