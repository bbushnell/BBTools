#!/bin/bash
# Read-only provenance audit: in=candidates.tsv labels=labels.tsv out=audit.tsv
# Stop codes are historical walk reasons; endpoint junction flags are later
# count queries and must not substitute for those reasons. No labels are changed.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKJoinAudit "$@"
