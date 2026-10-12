#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
for test_class in Prok16sReleaseTest Prok16sPartialRescueTest RrnaScoreOnlyTableTest;do
  java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" "prok.$test_class"
done
