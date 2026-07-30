#!/bin/bash
#Written by Noire
#Dev tool (not for normal users). Runs from releases/bbmap/dev/ -- finds classes in ../current/.

usage(){
echo "
Written by Noire

Description: A/B benchmark + correctness check for Clade.parseClade (positional) vs
             Clade.parseCladeFlex (order-independent, BytePrefixDispatcher). Verifies
             byte-identical Clades (canonical toBytes) and compares single-thread + full
             cold-load + multithreaded parse speed on a real spectra file.

Usage: cladeparsebench.sh in=<spectra.gz> [passes=8] [threads=N]

Java Parameters:
-Xmx            Set Java's memory usage, overriding autodetection (e.g. -Xmx24g).
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then usage; exit; fi

#Resolve the BBTools root (parent of dev/), following symlinks, so ../current/ is found.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")/.."
DIR="$(pwd)/"
popd > /dev/null

CP="$DIR""current/"

setEnv(){
  . "$DIR""javasetup.sh"
  . "$DIR""memdetect.sh"
  parseJavaArgs "--xmx=16g" "--xms=16g" "--percent=42" "--mode=auto" "$@"
  setEnvironment
}
setEnv "$@"

launch(){
  local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP clade.CladeParseBench $@"
  echo "$CMD" >&2
  eval $CMD
}
launch "$@"
