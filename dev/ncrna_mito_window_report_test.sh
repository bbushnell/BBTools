#!/usr/bin/env bash
# outdir=<fresh-directory> sha80=<canonical-executable>
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx1g -cp "$root/current" test.prok.NcrnaMitoWindowReportTest "$@"
