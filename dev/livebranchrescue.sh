#!/bin/bash
# Experimental single-K raw-contig comparison. Default: no network, native decisions.
# in=reads.fq.gz out=NEW.fa k=95 [net=model.bbnet cutoff=0.52372265]
# Forces raw-unitig and scalar inference settings; no reference/label input or production default change.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchRescue "$@"
