#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Yelan
Last modified October 8, 2026

Description:  Grows .bbnet neural networks to wider layer dimensions and
optionally deletes low-magnitude active edges.
Modified sparse output layers support inference and training with their stored input indices.

Usage:  modifynn.sh in=<old.bbnet> out=<grown.bbnet> dims=N,H1,H2,O

in=<file>       Input .bbnet network.
out=<file>      Output .bbnet network.
dims=<list>     Comma-separated absolute layer widths.  Omit for no-op/prune.
seed=1          Random seed for new active edges.
newweight=1e-3  Maximum absolute new-edge magnitude.
pruneabs=0      Delete old active edges with abs(weight) below this threshold.
zero2epsilon=f  Replace stored zero edges with seeded weights in the newweight range.
                Applies to legacy dense active zeros or explicit sparse zero entries.
                Dense zero-absent slots, absent sparse edges, and biases are unchanged.
privateperhead=0  Append this many private last-hidden nodes per new output head.
                Requires partition metadata and matching dims. With 0, new hidden nodes are shared.
report=<file>   Optional TSV modification report.
overwrite=f     (ow) Permit overwriting the output file.

Java Parameters:
-Xmx            Set Java memory usage; e.g. -Xmx1g.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter problems.
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2000m" "--xms=2000m" "--percent=10" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.ModifyNN $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" ml.ModifyNN "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
