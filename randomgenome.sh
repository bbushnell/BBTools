#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 11, 2025

Description:  Generates a random, (probably) repeat-free genome.

Usage:  randomgenome.sh len=<total size> chroms=<int> gc=<float> out=<file>

Parameters:
out=<file>      Output.
in=<file>       Optional input clade or fasta file.  If specified, the
                synthetic genome will conserve the input kmer frequencies.
k=5             Kmer length for base frequencies (2-5).
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
len=100000      Total genome size.
chroms=1        Number of pieces.
gc=0.5          GC fraction.
nopoly=f        Ban homopolymers.
pad=0           Add this many Ns to contig ends; does not count toward
                genome size.
seed=-1         Set to a positive number for deterministic output.
amino=f         Produce random amino acids instead of nucleotides.
includestop=f   Include stop codons in random amino sequences.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP synth.RandomGenome $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
