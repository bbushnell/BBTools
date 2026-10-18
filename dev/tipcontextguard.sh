#!/bin/bash
# EXPERIMENTAL: guard substitutions at one contig end; no assembly defaults changed.
# Trusted ConsensusMaker .alm files only, one graph each. originalmodel must use
# original unpadded alignments, NOT context-pooled support. Context must preserve
# the original target at zero-based offset. Assembly must include that named target.
# originalmodel=original.alm contextmodel=context.alm assembly=assembly.fa
# out=guarded.fa report=decisions.tsv offset=0 side=right minoverlap=60 maxoverlap=200
# Fixed experiment: mafsub=.5, mindepth=2, original allele supported twice is kept.
# One compatible end overlap within bounds, <=1 mismatch with 16bp anchors.
# Not proof of genomic uniqueness; a true one-read allele can still be overwritten.
# originalkmers=FILE optionally adds a mapping-independent veto from exact native
# KmerWindowDepth k=32 queries. Any original word with depth>=2 preserves the base.
# All overlapping original32mer queries must be present to allow a weak-site change.
# Missing evidence abstains; counts never force an edit. Default omits this veto.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --xss=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" consensus.TipContextGuard "$@"
