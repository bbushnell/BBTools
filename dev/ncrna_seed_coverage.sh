#!/usr/bin/env bash
# keysupport=t additionally writes .keys.tsv and .key_summary.tsv: distinct
# input-member support per loaded forward17-mer. Use fold=ALL representatives=f
# for an all-admitted-member audit; the default coverage outputs are unchanged.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.NcrnaSeedCoverage "$@"
