#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=3 -Xmx1g -cp "$root/current" test.prok.NcrnaAnnotatedCorpusTest "$@"
