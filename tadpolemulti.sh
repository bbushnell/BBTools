#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 9, 2026

Description:  Assemble at one kmer length, then process the ordered K list.
Each phase fuses eligible shorter-K overlaps and bridges remaining gaps,
sharing one count table.  Bridge K may be longer or shorter than assembly K.

Usage:  tadpolemulti.sh in=<reads> out=<contigs> k=96,124,64,32
Custom: tadpolemulti.sh in=<reads> out=<contigs> assemblek=96 fusek=64 bridgek=128,96,64,32 graphk=96

Core parameters:
k=             Ordered K list. The first value assembles by default; repeated
               values request repeated phases, e.g.96,124,64,32,96.
korder=input   input preserves order; legacy sorts unique Ks and runs all
               fusions before all bridges, with the old longest-K defaults.
assemblek=auto Initial assembly K.  An explicit value overrides k shorthand.
               The first matching list occurrence is consumed by assembly;
               if absent, assembly precedes the list.
fusek=auto     Exact reciprocal tip-overlap K values; all must be below assemblek.
               joink is an alias.  Set to none to disable.
bridgek=auto   Read-supported unbranched-walk K values; these may be above,
               below, or equal to assemblek.  Set to none to disable.
               In input order, fusek/bridgek filter the k list, which must
               contain all requested values. Without k, use fusek order then
               additional bridgek values after assembly.
graphk=auto    Final graph K for simplification, graph output or path extraction.
               Defaults to the last listed K; a matching final table is reused.
               A different graphk needs another load. Legacy defaults to assemblek.
fusenet=null   Optional fusion-join network; disabled by default. Supply the path
               to networks/tadpole_fusion.bbnet for the bundled dense model.
fusencutoff=   Required with fusenet. The bundled model was tested at 0.667098
               on simulated bacterial reads with minprob=0 minprobmain=f.
fusecoverageratio=0
               Maximum higher/lower whole-contig mean-depth ratio for exact
               overlap fusion, including final graph-K fusions. Not gap bridges.
               0 disables; otherwise finite >=1. Zero-depth ends then reject.
graphmergecoverageratio=0
               Same ratio check for ordinary direct merges in the final graph.
               Requests final graph processing, reusing a matching last table.
               Does not filter initial assembly, cross-K fusion/bridging,
               indirect bubble removal, or path extraction. 0 disables;
               otherwise finite >=1. Zero-depth ends then reject.
               Both guards work without an NN; 1.75 was tested with the bundled
               NN. They can reject correct joins across genuine depth changes.
crosskmaxdepthratio=3
               Reject a bridge whose depth exceeds this multiple of the
               greater flank coverage.  Set to 0 to disable.
crosskpasses=10
               Maximum direct-merge passes after each fuse or bridge phase.
crosskmaxlen=500
               Maximum unbranched bridge distance searched from an eligible
               contig end.  Cycles and branches terminate earlier.

Optional coverage safeguards: fusecoverageratio=1.75 graphmergecoverageratio=1.75
These use stored contig depths, not an additional evidence table.

Other Tadpole parameters, including pop, shave, rinse, mincountseed, and
mincountextend, are passed to the initial assembly.
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then usage; exit; fi

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

CP="$DIR""current/"

setEnv(){
  . "$DIR""javasetup.sh"
  . "$DIR""memdetect.sh"
  parseJavaArgs "--xmx=14g" "--xms=14g" "--percent=84" "--mode=auto" "$@"
  setEnvironment
}
setEnv "$@"

launch(){
  local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP assemble.TadpoleMulti $@"
  echo "$CMD" >&2
  java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.TadpoleMulti "$@"
}
launch "$@"
