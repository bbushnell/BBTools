#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx4g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=3 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.NcrnaMitoSsuCompareTest "$@"
