#!/usr/bin/env bash
# Run source-coordinate and tile-conservation fixtures in a fresh directory.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" test.prok.NcrnaMitoContigMappingSummaryTest "$@"
