#!/usr/bin/env bash
# Native prok truth/grade fixture and6S subtype/union controls.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx1g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.NcrnaFreshPanelGradeTest "$@"
