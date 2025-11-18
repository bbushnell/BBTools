#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 27, 2018

Description:  Generates a prokaryotic gene model (.pkm) for gene calling.
Input is fasta and gff files.
The .pkm file may be used by CallGenes.

Usage:  analyzegenes.sh in=x.fa gff=x.gff out=x.pgm

File parameters:
in=<file>       A fasta file or comma-delimited list of fasta files.
gff=<file>      A gff file or comma-delimited list.  This is optional;
                if present, it must match the number of fasta files.
                If absent, a fasta file 'foo.fasta' will imply the
                presence of 'foo.gff'.
out=<file>      Output pgm file.

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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.AnalyzeGenes $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
