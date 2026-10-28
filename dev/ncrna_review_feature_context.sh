#!/usr/bin/env bash
# Exact original-GFF context for prediction-derived loci; no CM or truth decisions.
# Producer verifies gffhashes and fastahashes against source bytes before invocation.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" prok.NcrnaReviewFeatureContext "$@"
