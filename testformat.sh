#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 6, 2025

Description:  Tests file extensions and contents to determine format,
quality, compression, interleaving, and read length.  More than one file
may be specified.  Note that ASCII-33 and ASCII-64 cannot always
be differentiated.

Usage:  testformat.sh <file>

See also:  testformat2.sh, stats.sh

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

	parseJavaArgs "--xmx=120m" "--xms=120m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP fileIO.FileFormat $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
