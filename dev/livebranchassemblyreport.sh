#!/bin/bash
# Summarizes a collected native/live VALID panel, preserving empty/undefined cases.
# in=directory out=NEW_directory expected=153
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchAssemblyReport "$@"
