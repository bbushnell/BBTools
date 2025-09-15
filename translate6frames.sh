#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Translates nucleotide sequences to all 6 amino acid frames,
or amino acids to a canonical nucleotide representation.
Input may be fasta or fastq, compressed or uncompressed.

Usage:  translate6frames.sh in=<input file> out=<output file>

See also:  callgenes.sh

Optional parameters (and their defaults)

Input parameters:
in=<file>       Main input. in=stdin.fa will pipe from stdin.
in2=<file>      Input for 2nd read of pairs in a different file.
int=auto        (interleaved) t/f overrides interleaved autodetection.
qin=auto        Input quality offset: 33 (Sanger), 64, or auto.
aain=f          False if input is nucleotides, true for amino acids.
reads=-1        If positive, quit after processing X reads or pairs.

Output parameters:
out=<file>      Write output here.  'out=stdout.fa' goes to standard out.
out2=<file>     Use this to write 2nd read of pairs to a different file.
overwrite=t     (ow) Grant permission to overwrite files.
append=f        Append to existing files.
ziplevel=2      (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80    Length of lines in fasta output.
qout=auto       Output quality offset: 33 (Sanger), 64, or auto.
aaout=t         False to output nucleotides, true for amino acids.
tag=t           Tag read id with the frame, adding e.g. ' fr1'
frames=6        Only print this many frames.  
                If you already know the sense, set 'frames=3'.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

z="-Xmx2g"
z2="-Xms2g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 2000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

translate6frames() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.TranslateSixFrames $@"
	echo $CMD >&2
	eval $CMD
}

translate6frames "$@"
