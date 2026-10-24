#!/usr/bin/env bash
# Export source-bound call neighborhoods; all coordinates stay forward-source.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=4 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" prok.NcrnaMitoCallContexts "$@"
