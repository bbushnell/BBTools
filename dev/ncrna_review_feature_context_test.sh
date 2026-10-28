#!/usr/bin/env bash
# Arguments: fresh fixture directory, canonical sha80 executable.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx1g -cp "$root/current" test.prok.NcrnaReviewFeatureContextTest "$@"
