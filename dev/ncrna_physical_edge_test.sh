#!/usr/bin/env bash
# Mechanism fixtures, not a calibration or a default change.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
for suite in NcrnaPhysicalEdgeTest Prok16sPartialRescueTest RrnaItsDeriverTest;do
  java -ea -Xmx2g -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" "test.prok.$suite" "$@"
done
