#!/bin/bash

#Running callvariants2.sh is equivalent to running callvariants.sh with the "multi" flag.
#See callvariants.sh for usage information.

#callvariants2 is intended for multiple sam/bam files, one from each sample, which should have variants called independently; the point is that allele frequencies will be reported for ALL samples at locations where ANY sample has a variant called.
#callvariants2 is NOT a better version of callvariants, it's the same, just designed for multisample processing.
#If you have only 1 sample (regardless of how many sam/bam files there are) you should use callvariants.sh without the "multi" flag.

usage(){
bash "$DIR"callvariants.sh
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.CallVariants2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
