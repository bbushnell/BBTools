#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx512m -cp "$root/current" test.prok.SeedModelPositionTableTest "$@"
