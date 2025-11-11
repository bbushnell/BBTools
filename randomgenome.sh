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

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

randomgenome() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP synth.RandomGenome $@"
	echo $CMD >&2
	eval $CMD
}

randomgenome "$@"
