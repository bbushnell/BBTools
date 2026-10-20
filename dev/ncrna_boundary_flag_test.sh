#!/usr/bin/env bash
# Fresh-JVM CLI-to-loaded-boundary fixture; pass expectations then CallGenes flags.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" test.prok.NcrnaBoundaryFlagTest "$@"
