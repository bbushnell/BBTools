#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
for suite in GeneCallerSsuScoreTest CallGenesRrnaFallbackConfigTest;do
  java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx1g -cp "$root/current" "test.prok.$suite"
done
