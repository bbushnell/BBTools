#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 10, 2026

Description:  Evaluates a MAG-QC 'denominator' subnet by the quality of the
DERIVED subset-relative completeness, not just its raw regression target.
Such a subnet predicts a gene subset's native COMPLEMENT (e.g. an organism's
total ncRNA count); the useful quantity is completeness = observed /
predicted-native.  Reports the raw native-count R2 (how well the denominator
itself is predicted) AND the derived-completeness R2/MAE/RMSE
(observed/predicted-native vs observed/true-native, clamped to [0,1]).

This is a DEVELOPMENT tool (subnet evaluation); it is not part of the
end-user MAG-QC path.

Usage:  subnetratioscore.sh net=<subnet.bbnet> in=<val_vectors.tsv> [numobs=5]

Required parameters:
net=<file>      Trained subnet .bbnet (RegressionTrainer output).
in=<file>       Subnet vector TSV with '#dims' header; last column is the TRUE
                native complement.
Optional parameters (and their defaults):
numobs=5        Number of leading observed-count columns (their sum = the
                observed numerator).

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

	parseJavaArgs "--xmx=8g" "--xms=1g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.SubnetRatioScore $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
