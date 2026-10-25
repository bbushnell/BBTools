#!/usr/bin/env bash
# Native fixture; supply a new scratch directory on Dori.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx1g -cp "$root/current" test.prok.NcrnaMitoUnmatchedEvidenceTest "$@"
