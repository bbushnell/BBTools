#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Selects reads with designated numeric IDs.
Paired mates write to out2 when set, otherwise immediately after selected reads.

Usage:  getreads.sh in=<file> id=<number,number,number...> out=<file>

The first read (or pair) has ID 0, the second read (or pair) has ID 1, etc.

Parameters:
in=<file>       Specify the input file, or stdin.
in2=<file>      Optional second input file for paired reads.
qfin=<file>     Optional QUAL sidecar for nucleotide FASTA input.
qfin2=<file>    Optional QUAL sidecar for the second FASTA input.
                QUAL names are literal; # expansion is not applied.
                Requires nucleotide FASTA with minscaf=1, no fastareadlen
                splitting, and forcesectionname=f.
                Supports unpaired input, or separate paired in/in2 files.
out=<file>      Specify the output file, or stdout. # becomes 1 and 2.
out2=<file>     Optional output for selected mates; otherwise mates write to out.
qfout=<file>    Optional QUAL sidecar for selected reads in FASTA out.
qfout2=<file>   Optional QUAL sidecar for selected mates in FASTA out2.
                QUAL output names are literal; # expansion is not applied.
                stdout and stderr are not supported for QUAL sidecars.
id=             Comma delimited list of numbers or ranges, in any order.
                For example:  id=5,93,17-31,8,0,12-13

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.GetReads $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.GetReads "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"