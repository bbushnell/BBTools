#!/usr/bin/env bash
# Join frozen source-annotation evidence; no mapping or CM search.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=4 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" prok.NcrnaMitoSsuAnnotationJoin "$@"
