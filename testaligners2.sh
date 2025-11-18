#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 31, 2025

Description:  Tests multiple aligners using random sequences.
The sequences have variable pairwise ANI, and each
ANI level is tested multiple times for average accuracy
and loop count.
Outputs the identity, rstart and rstop positions, time, and #loops.
Note that the 'design' ANI is approximate and will not match
the measured ANI.

Usage:
testaligners2.sh iterations=30 maxani=100 minani=90 step=2

Parameters:
length=40k      Length of sequences.
iterations=32   Iterations to average; higher is more accurate.
maxani=80       Max ANI to model.
minani=30       Min ANI to model.
step=2          ANI step size.
sinewaves=0     Sinewave count to model variable conservation.
threads=        Parallel alignments; default is logical cores.
simd            Enable SIMD operations; requires AVX-256 and Java 17+.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$0"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.TestAlignerSuite $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
