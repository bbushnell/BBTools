#!/usr/bin/env bash
# Builds the two reference indexes used by experimental protein HBM assignment.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: hbmreferenceindex.sh tool=coveringset|sidecar [producer arguments] -Xmx8g
coveringset: families=TABLE out=FILE summary=FILE alphabet=amino k=5 kdesign=6
  step=50 target=.999 minhits=1 maxfamilies=0 streamfamilies=t t=1
sidecar: consensus=FASTA familiesdir=NUMERIC_FASTA_DIRECTORY kmersets=FILE
  out=FILE variant=centered rawpseudo=0 background=families
The bare argument selftest runs the selected producer's existing fixture.
This launcher does not activate the experimental reference in ProkCC.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=(); TOOL=
for arg in "$@"; do
 case ${arg,,} in
  tool=*) [[ -z $TOOL ]] || { echo 'Duplicate tool argument' >&2; exit 2; }; TOOL=${arg#*=};;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
case ${TOOL,,} in
 coveringset) CLASS=prok.CoveringSet;;
 sidecar|familyshortlistsidecarbuilder) CLASS=prot.FamilyShortlistSidecarBuilder;;
 *) echo "Unknown or missing tool: $TOOL" >&2; exit 2;;
esac
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
