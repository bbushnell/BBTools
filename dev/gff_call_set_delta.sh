#!/usr/bin/env bash
# Compare sorted GFF call sets: before=FILE after=FILE outprefix=PREFIX.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.GffCallSetDelta "$@"
