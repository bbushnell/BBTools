#!/usr/bin/env bash
# Test direct-loader/caller agreement using one explicit TPP fixture directory.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 \
  -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" test.prok.CallGenesTppTest "$@"
