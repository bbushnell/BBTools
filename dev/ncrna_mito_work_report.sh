#!/usr/bin/env bash
# Bound saved-log tier accounting only; no caller or alignment execution.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx1g -cp "$root/current" prok.NcrnaMitoWorkReport "$@"
