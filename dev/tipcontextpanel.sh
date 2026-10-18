#!/bin/bash
# Experimental bounded panel orchestration; not ordinary assembly/polishing.
# mode=prepare in=assembly.fa targets=selected.bed dir=NEW_DIRECTORY
# mode=merge in=assembly.fa dir=PANEL_DIRECTORY out=guarded.fa report=edits.tsv
# mode=grade in=edits.tsv before=original.used_snps after=guarded.used_snps out=summary.tsv
# prepare uses the frozen 60..200bp overlap rule and adds unique donor controls.
# merge requires every prepared side to have side.guarded.fa; conflicting edits abstain.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" consensus.TipContextPanel "$@"
