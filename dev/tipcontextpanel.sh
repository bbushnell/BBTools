#!/bin/bash
# Experimental bounded panel orchestration; not ordinary assembly/polishing.
# mode=prepare in=assembly.fa targets=selected.bed dir=NEW_DIRECTORY
# mode=merge in=assembly.fa dir=PANEL_DIRECTORY out=guarded.fa report=edits.tsv
# mode=grade in=edits.tsv before=original.used_snps after=guarded.used_snps out=summary.tsv
# prepare defaults to minoverlap=60 maxoverlap=200 and adds unique donor controls.
# Experimental minoverlap>=33 retains16bp mismatch anchors. Use identical bounds
# for the later tipcontextguard.sh call; these do not relax its evidence checks.
# indexed=t optionally uses an exact16mer prefix index; default=f retains exhaustive discovery.
# merge requires every prepared side to have side.guarded.fa; conflicting edits abstain.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" consensus.TipContextPanel "$@"
