#!/usr/bin/env bash
# Saved-call grading against explicitly reviewed fresh prok CM truth.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx4g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaFreshPanelGrade "$@"
