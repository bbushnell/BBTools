#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2019

Description:  Writes a shell script to download one genome assembly
and gff per genus or species, from ncbi.  Attempts to select the best assembly
on the basis of contiguity.

Usage:  fetchproks.sh <url> <outfile> <max species per genus: int> <use best: t/f>

Examples:
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/ bacteria.sh 2 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/archaea/ archaea.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/ viral.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/protozoa/ protozoa.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/invertebrate/ invertebrate.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/fungi/ fungi.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/plant/ plant.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/ vertebrate_mammalian.sh 0 true
fetchproks.sh ftp://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_other/ vertebrate_other.sh 0 true

Mitochondrion, plasmid, and plastid are different and use gbff2gff.

Processing parameters:
None yet!

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

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.FetchProks $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
