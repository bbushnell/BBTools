#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 6, 2024

Description:  Converts BGI/Complete Genomics reads to Illumina header format,
and optionally appends barcodes/indexes. For example, 
@E200008112L1C001R00100063962/1 
would become
@E200008112:0:FC:1:6396:1:1 1:N:0:

Usage:  cg2illumina.sh in=<input file> out=<output file> barcode=<string>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary output, or read 1 output.
out2=<file>     Read 2 output if reads are in two files.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
barcode=        (index) Optionally append a barcode to the header.
parseextra=f    Set this to true if the reads headers have comments 
                delimited by a whitespace.

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

z="-Xmx300m"
z2="-Xms300m"
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

function bgi2ill() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP hiseq.BGI2Illumina $@"
	echo $CMD >&2
	eval $CMD
}

bgi2ill "$@"
