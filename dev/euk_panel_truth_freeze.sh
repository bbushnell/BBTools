#!/usr/bin/env bash
# Invoke only with the review/authorization bindings checked by the producer.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx1g -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" prok.EukPanelTruthFreeze "$@"
