#!/usr/bin/env bash
# Compare modern/legacy18S resource bindings against a retained installation.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" prok.Euk18sFamilyResourceTest "$@"
