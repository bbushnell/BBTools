#!/usr/bin/env bash
# outdir=<fresh-directory> sha80=<canonical-executable>
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=4 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" test.prok.NcrnaMitoCoverageReportTest "$@"
