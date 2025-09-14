#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description: Filters barcodes by quality, and generates quality histograms.

Usage:       filterbarcodes.sh in=<file> out=<file> maq=<integer>

Input parameters:
in=<file>       Reads that have already been muxed with barcode qualities using mergebarcodes.sh.
int=auto        (interleaved) If true, forces fastq input to be paired and interleaved.
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Output parameters:
out=<file>      Write filtered reads here.  'out=stdout.fq' will pipe to standard out.
cor=<file>      Correlation between read and index qualities.
bqhist=<file>   Barcode quality histogram by position.
baqhist=<file>  Barcode average quality histogram.
bmqhist=<file>  Barcode min quality histogram.
overwrite=t     (ow) Set to false to force the program to abort rather than overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
fastawrap=80    Length of lines in fasta output.
qout=auto       ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
maq=0           Filter reads with barcode average quality less than this.
mmq=0           Filter reads with barcode minimum quality less than this.

Other parameters:
pigz=t          Use pigz to compress.  If argument is a number, that will set the number of pigz threads.
unpigz=t        Use pigz to decompress.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

filterbarcodes() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.CorrelateBarcodes $@"
	echo $CMD >&2
	eval $CMD
}

filterbarcodes "$@"
