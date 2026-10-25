#!/usr/bin/env bash
# Saved calls, inclusive source203 and provisional raw CM hits; no truth freeze.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx8g -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaClippingComparison "$@"
