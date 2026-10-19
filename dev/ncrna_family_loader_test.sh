#!/usr/bin/env bash
# Run migration fixtures against baseline=/absolute/pre-migration/resources.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" test.prok.NcrnaFamilyLoaderTest "$@"
