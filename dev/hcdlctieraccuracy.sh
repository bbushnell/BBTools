#!/bin/bash
#Written by Chloe
#April 14, 2026

usage(){
echo "
Written by Chloe
Last modified April 14, 2026

Description:  Measures per-tier DLC and HC accuracy as a function of
              V (DLC), hcBeff (HC), and hcUnseen (HC) for UDLL6.  Emits
              three TSV files: {prefix}dlc_v.tsv, {prefix}hc_beff.tsv,
              {prefix}hc_unseen.tsv.  Each file has linAvg, geoAvg, and
              harmAvg absolute-error columns so all three averaging modes
              are produced in a single run.

Usage:  hcdlctieraccuracy.sh hbits=2 ddls=100000 maxmult=8192

Flags:
  hbits=2        History bits.  0 = DLC only; 1 or 2 use UDLL6; 3 reserved.
  buckets=2048   Bucket count B (fixed experiment).
  ddls=100000    Number of independent DDLs to average.
  maxmult=8192   Maximum cardinality in units of B.
  sample=logcard Sampling mode: logcard | entry | entryexit | add.
  points=300     Log-spaced checkpoints (logcard mode only).
  out_prefix=hcdlc_  Output file prefix.
  threads=auto   Worker thread count; one DDL per thread at a time.
  seed=-1        Master seed (-1 = use nanoTime).

Java Parameters:
-Xmx          This will set Java's memory usage, overriding autodetection.
              -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
              specify 200 megs.  The max is typically 85% of physical memory.
-eoom         This flag will cause the process to exit if an out-of-memory
              exception occurs.  Requires Java 8u92+.
-da           Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")/.."
DIR="$(pwd)/"
popd > /dev/null

#DIR="$(dirname "$(readlink -f "$0")")/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

hcdlc() {
	local CMD="java $EA $EOOM $z -cp $CP cardinality.HCDLCTierAccuracy $@"
	echo $CMD >&2
	eval $CMD
}

hcdlc "$@"
