#!/bin/bash
# Read-only exact32bp origin audit; reference never feeds the joining decision.
# in=assembly.fa ref=reference.fa positions=contig_index:start0,... out=hits.tsv
# At most32 selected arms; all coordinates/indices are zero-based. test=t runs tiny checks.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowTipAnchorAudit "$@"
