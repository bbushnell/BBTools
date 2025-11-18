#!/bin/bash

usage(){
echo "
bbmerge-auto.sh is a wrapper for BBMerge that attempts to use all available
memory, instead of a fixed amount.  This is for use with the Tadpole options
of error-correction (ecct) and extension, which require more memory.
For merging by overlap only, please use bbmerge.sh.  If you set memory
manually with the -Xmx flag, bbmerge.sh and bbmerge-auto.sh are equivalent.

For information about usage and parameters, please run bbmerge.sh.
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

	parseJavaArgs "--xmx=15000m" "--xms=15000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.BBMerge $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
