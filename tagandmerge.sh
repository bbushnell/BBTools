#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 16, 2024

Description:  Accepts multiple input files from a demultiplexed lane.
Parses the barcode from the filename and adds (tab)BARCODE to read headers.
Outputs all reads into a single file.  Optionally, trims bases and drops R2.
Intended for evaluating demultiplexing methods.  For example:
tagandmerge.sh path/*0.*.fastq.gz dropr2 trim out=tagged.fq.gz barcodes=bc.txt

Usage:  tagandmerge.sh *.fastq.gz out=<output file>
or
tagandmerge.sh in=<file,file,file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file,file>  A comma-delimited list of files.  If wildcards are used,
                omit in= and the commas.
out=<file>      Print all reads to this destination.
barcodes=<file> Print barcodes from file names to this destination.
trim=-1         If positive, trim all reads to this length.
dropr2=f        Discard read 2 if the input is interleaved.
shrinkheader=f  (shrink) Illumina only; remove unnecessary header fields.
remap=-+        Remap symbols in the barcode.  By default, '+' replaces '-'.
                To eliminate this set 'remap=null'.

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

function tagandmerge() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP barcode.TagAndMerge $@"
	echo $CMD >&2
	eval $CMD
}

tagandmerge "$@"
