#!/bin/bash
#Written by Neptune
#Dev tool (not for normal users). Runs from releases/bbmap/dev/ -- finds classes in ../current/.

usage(){
echo "
Written by Neptune

Description: Enriches a ranking vector file by duplicating queries where the net's prediction fails (hard-example oversampling).

Usage: enrichfailures.sh <args>   (arguments are parsed by ml.EnrichFailures; run the class or read
        its source for the exact flags -- this is a developer/analysis tool.)

Java Parameters:
-Xmx            Set Java's memory usage, overriding autodetection (e.g. -Xmx8g).
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
  parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
  setEnvironment
}
setEnv "$@"

launch(){
  local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.EnrichFailures $@"
  echo "$CMD" >&2
  eval $CMD
}
launch "$@"
