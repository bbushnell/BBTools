#!/usr/bin/env bash
# Explicit pilot/full/select phases; frozen source bindings required.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=16 -XX:+ExitOnOutOfMemoryError -Xmx8g -cp "$root/current" prok.NcrnaMitoCoverageSelect "$@"
