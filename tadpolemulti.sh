#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 30, 2026

Description:  Assemble at one kmer length, directly join unique reciprocal
tip overlaps at requested shorter K values, then load independent bridge-K
tables and conservatively bridge remaining gaps.  Bridge K may be longer or
shorter than the assembly K.

Usage:  tadpolemulti.sh in=<reads> out=<contigs> k=31,75,140
Custom: tadpolemulti.sh in=<reads> out=<contigs> assemblek=96 fusek=64 bridgek=128,96,64,32 graphk=96

Core parameters:
k=             Shorthand kmer list; internally sorted descending.  The longest
               value assembles, and shorter values fuse and bridge by default.
assemblek=auto Initial assembly K.  An explicit value overrides k shorthand.
fusek=auto     Exact reciprocal tip-overlap K values; all must be below assemblek.
               joink is an alias.  Set to none to disable.
bridgek=auto   Read-supported unbranched-walk K values; these may be above,
               below, or equal to assemblek.  Set to none to disable.
graphk=auto    Final graph K for simpleomnitigs or graphcover.  Defaults to
               assemblek; a matching final bridge table is reused when available.
crosskmaxdepthratio=3
               Reject a bridge whose depth exceeds this multiple of the
               greater flank coverage.  Set to 0 to disable.
crosskpasses=10
               Maximum direct-merge passes after each fuse or bridge phase.
crosskmaxlen=500
               Maximum unbranched bridge distance searched from an eligible
               contig end.  Cycles and branches terminate earlier.

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
  eval $CMD
}
launch "$@"
