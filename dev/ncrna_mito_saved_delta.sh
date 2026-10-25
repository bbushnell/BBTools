#!/usr/bin/env bash
# Saved audit/reference and historical-evidence joins only.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.NcrnaMitoSavedDelta "$@"
