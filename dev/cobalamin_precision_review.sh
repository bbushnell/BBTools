#!/usr/bin/env bash
# Joins frozen CM evidence, explicit decisions and saved calls; never changes truth.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.CobalaminPrecisionReview "$@"
