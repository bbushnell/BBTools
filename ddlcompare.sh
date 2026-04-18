#!/bin/bash

usage(){
echo "
Written by Ady
Last modified April 17, 2026

Description:  Compares two genomes using DynamicDemiLog bucket matching
to estimate WKID and ANI from k-mer intersection.

Usage:  ddlcompare.sh <genome1.fa> <genome2.fa> [k=31] [buckets=2048]

Parameters:
k=31            K-mer length for hashing.
buckets=2048    Number of DDL buckets.

Example:
ddlcompare.sh ecoli.fa mruber.fa
ddlcompare.sh ref.fa mutant.fa k=31 buckets=2048
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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLCompare $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
