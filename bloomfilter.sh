#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 20, 2022

Description:  Filters reads potentially sharing a kmer with a reference.
The more memory, the higher the accuracy.  Reads going to outu are guaranteed
to not match the reference, but reads going to outm might may or may not
match the reference.

Usage:  bloomfilter.sh in=<input file> out=<nonmatches> outm=<matches> ref=<reference>

Example:
bloomfilter.sh in=reads.fq outm=nonhuman.fq out=human.fq k=31 minhits=3 ref=human.fa

Error correction and depth filtering can be done simultaneously.

File parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
outm=<file>     (out) Primary matched read output.
outm2=<file>    (out2) Matched read 2 output if reads are in two files.
outu=<file>     Primary unmatched read output.
outu2=<file>    Unmatched read 2 output if reads are in two files.
outc=<file>     Optional output stream for kmer counts.
ref=<file>      Reference sequence file, or a comma-delimited list.
                For depth-based filtering, set this to the same as the input.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Hashing parameters:
k=31            Kmer length.
hashes=2        Number of hashes per kmer.  Higher generally reduces 
                false positives at the expense of speed.
sw=t            (symmetricwrite) Increases accuracy when bits>1 and hashes>1.
minprob=0.5     Ignore reference kmers with probability of being correct
                below this (affects fastq references only).
memmult=1.0     Fraction of free memory to use for Bloom filter.  1.0 should
                generally work; if the program crashes with an out of memory
                error, set this lower.  Higher increases specificity.
cells=          Option to set the number of cells manually.  By default this
                will be autoset to use all available memory.  The only reason
                to set this is to ensure deterministic output.
seed=0          This will change the hash function used.
bits=           Bits per cell; it is set automatically from mincount.

Reference-matching parameters:
minhits=3       Consecutive kmer hits for a read to be considered matched.
                Higher reduces false positives at the expense of sensitivity.
mincount=1      Minimum number of times a read kmer must occur in the 
                reference to be considered a match (or printed to outc).
requireboth=f   Require both reads in a pair to match the ref in order to go
                to outm.  By default, pairs go to outm if either matches.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bloomfilter() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP bloom.BloomFilterWrapper $@"
	echo $CMD >&2
	eval $CMD
}

bloomfilter "$@"
