#!/bin/bash
# Summarizes a collected assembly panel, preserving empty/undefined cases.
# in=directory out=NEW_directory expected=153 fold=validation
# fold=confirmation explicitly requires fresh confirmation labels on every row.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchAssemblyReport "$@"
