#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
# Worker count is explicit (t=, default1); keep GC/JIT overhead bounded too.
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=16 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.NcrnaResourceReport "$@"
