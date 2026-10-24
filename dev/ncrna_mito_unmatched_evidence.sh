#!/usr/bin/env bash
# Saved native audit/GFF/reference associations; full contigs stay on the cluster.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx4g -cp "$root/current" prok.NcrnaMitoUnmatchedEvidence "$@"
