#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified August 21, 2026

Description:  Converts a TaxTree between portable TSV and legacy serialized
formats according to the output extension, reloads it, and validates every
persisted node and tree field.

Usage:  taxtree2.sh in=tree.taxtree.gz out=taxtree.tsv.gz

Parameters:
in=<file>       TaxTree input.  Defaults to the bundled tree.
out=<file>      Output; .tsv[.gz] selects text, otherwise legacy serialization.
overwrite=f     Overwrite an existing output file.

Java Parameters:
-Xmx            Set Java memory usage, overriding autodetection.
-eoom           Exit on an out-of-memory exception.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
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

launch(){
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP tax.TaxTreeText $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
