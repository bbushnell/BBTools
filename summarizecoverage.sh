#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 5, 2020

Description:  Summarizes coverage information from basecov files
              created by pileup.sh.  They should be named like
              'sample1_basecov.txt' but other naming styles are fine too.

Usage:        summarizecoverage.sh *basecov.txt out=<output file>

Parameters:
in=<file>           'in=' is not necessary.  Any filename used as a
                    parameter will be assumed to be an input basecov file.
out=<file>          Write the summary here.  Default is stdout.
reflen=-1           If positive, use this as the total reference length.
                    Otherwise, assume basecov files report every ref base.

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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP covid.SummarizeCoverage $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
