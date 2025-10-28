#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 14, 2025

Description:  Counts the number of unique kmers in a file.
Prints a fasta or tsv file containing all kmers and their counts.
Supports K=1 to 15, though values above 8 should use KmerCountExact.
SEE ALSO: kmercountexact.sh

Usage:   kmercountshort.sh in=<file> out=<file> k=4

Input may be fasta or fastq, compressed or uncompressed.
Output may be stdout or a file.  out, khist, and peaks are optional.


Input parameters:
in=<file>           Primary input file.
in2=<file>          Second input file for paired reads.

Output parameters:
out=<file>          Print kmers and their counts.  Extension sensitive;
                    .fa or .fasta will produce fasta, otherwise tsv.
mincount=0          Only print kmers with at least this depth.
reads=-1            Only process this number of reads, then quit (-1 means all).
rcomp=t             Store and count each kmer together and its reverse-complement.
comment=            Denotes start of the tsv header.  E.g. 'comment=#'
skip=1              Count every Nth kmer.  If skip=2, count every 2nd kmer, etc.

Counting parameters:
k=4                 Kmer length - needs at least (threads+1)*8*4^k memory.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.
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

z="-Xmx1g"
z2="-Xms1g"
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

kmercountshort() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP jgi.KmerCountShort $@"
	echo $CMD >&2
	eval $CMD
}

kmercountshort "$@"
