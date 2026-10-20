#!/usr/bin/env bash
# Per-subtype 6S boundary override fixture; no files or caller resources required.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx1g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.SixsBoundaryConfigTest "$@"
