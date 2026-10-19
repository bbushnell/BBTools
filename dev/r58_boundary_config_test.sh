#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx1g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.R58BoundaryConfigTest "$@"
