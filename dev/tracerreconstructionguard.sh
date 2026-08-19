#!/bin/bash
#Written by Neptune
#Dev tool (not for normal users). Runs from releases/bbmap/dev/ -- finds classes in ../current/.

usage(){
echo "
Written by Neptune

Description: Regression guard for the idaligner Tracer traceback-reconstruction bug
(fixed Aug 19 2026). Sweeps all 7 idaligner aligners with a traceback path across a
substitution/indel case set and asserts the reconstructed matchString is
byte-for-byte consistent with the real query/ref sequences. Takes no arguments.

Usage: tracerreconstructionguard.sh

Exits 0 and prints ALL PASS if every case passes for every aligner.
Exits 1 with FAIL detail on stderr if any case fails -- treat as a real regression.

Java Parameters:
-Xmx            Set Java's memory usage, overriding autodetection (e.g. -Xmx1g).
-da             Disable assertions.
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then usage; exit; fi

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
  parseJavaArgs "--xmx=1g" "--xms=1g" "--percent=42" "--mode=auto" "$@"
  setEnvironment
}
setEnv "$@"

launch(){
  local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP idaligner.TracerReconstructionGuard $@"
  echo "$CMD" >&2
  eval $CMD
}
launch "$@"
