#!/usr/bin/env bash
# Read-only census of representative-grouped MMseqs membership TSVs.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: clustersizehistogram.sh in=clusters.tsv reps=representatives.faa out=PREFIX
  expectedfamilies=N expectedgenes=N optionally enforce independently known totals.
Input: exactly representative<TAB>member per row, grouped by representative.
Every FASTA representative must have one contiguous group containing itself once.
Coverage counts membership rows, not unique sequences or distinct organisms.
Outputs: PREFIX.histogram.tsv, .bins.tsv, .targets.tsv, .top.tsv, .summary.tsv.
Targets use the minimum family count within tied size bins. Existing outputs fail.
Runs one input-processing thread. Default heap is 4g; override with -Xmx.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case ${arg,,} in
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.ClusterSizeHistogram "${ARGS[@]}"
