#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx512m -XX:ActiveProcessorCount=2 --add-modules jdk.incubator.vector -cp "$root/current" test.clade.CladeBuilderTest "$@"
