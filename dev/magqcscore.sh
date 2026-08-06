#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 6, 2026

Description:  Scores a trained MAG-QC (CheckM3) regression net (.bbnet) on a
held-out validation TSV.  For each output it reports MAE and RMSE, the
predict-the-mean baseline MAE/RMSE, and R2 = 1 - net_SSE/baseline_SSE.  The
validation TSV must have the same '#dims' layout the net was trained on (as
written by magqcvectormaker.sh).

This is a DEVELOPMENT tool (model evaluation); it is not part of the end-user
MAG-QC path.

Note:  the output rows are labeled 'completeness' and 'contamination'; a
single-output net's one row is labeled 'completeness' regardless, so for a
contamination-only net that value IS the contamination metric.

Usage:  magqcscore.sh net=<net.bbnet> in=<val.tsv>

Required parameters:
net=<file>      Trained network in BBNet format.
in=<file>       Held-out validation vector TSV (with a '#dims' header).

Java Parameters:
-Xmx            Set Java heap (overrides autodetection).
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=1g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.MagQCScore $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
