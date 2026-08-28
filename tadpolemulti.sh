#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 23, 2026

Description:  Assemble at a long kmer length, directly join unique reciprocal
overlaps between low-depth contig ends at every requested shorter k, then load
the short-k tables and conservatively bridge remaining gaps.

Usage:  tadpolemulti.sh in=<reads> out=<contigs> k=31,75,140

Core parameters:
k=             Comma-delimited kmer lengths; internally sorted descending.
graphk=auto    Final graph K for simpleomnitigs or graphcover.  Defaults to
               the shortest K, whose table is reused when available.
crosskmaxdepthratio=3
               Reject a bridge whose depth exceeds this multiple of the
               greater flank coverage.  Set to 0 to disable.
crosskpasses=10
               Maximum direct-merge passes at each shorter kmer length.
crosskmaxlen=500
               Maximum unbranched short-k distance searched from a low-depth
               contig end.  Cycles and branches terminate earlier.

Other Tadpole parameters, including pop, shave, rinse, mincountseed, and
mincountextend, are passed to the initial long-k assembly.
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
