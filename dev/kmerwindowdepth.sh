#!/bin/bash
# Exact unfiltered read counts along forward-reference windows; one table per run.
# in=reads ref=fasta scaffold=record_name regions=first-last[,first-last] out=tsv k=N
# Region positions are ONE-BASED INCLUSIVE KMER ENDS, not windows of kmer starts.
# The named reference is query-only. It never contributes to the read counts.
# Outputs the central depth and A/C/G/T extensions on each side, in reference orientation.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.KmerWindowDepth "$@"
