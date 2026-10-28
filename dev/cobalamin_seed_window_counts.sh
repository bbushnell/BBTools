#!/usr/bin/env bash
# Saved producing-window distinct17-mer histogram; no caller or alignment replay.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.CobalaminSeedWindowCounts "$@"
