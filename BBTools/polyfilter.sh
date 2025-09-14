#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 1, 2024

Description:  Filters reads to remove those with suspicious homopolymers.

Usage:  polyfilter.sh in=<input reads> out=<filtered reads>

Example:
polyfilter.sh in=reads.fq out=clean.fq outb=bad.fq k=31 polymers=G


File parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Output for clean reads.
outb=<file>     Output for bad (homopolymer) reads.
extra=<file>    Comma-delimited list of additional sequence files.
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

Filtering rules:
Reads will always be discarded if they fails ldf2, entropy2, or minpolymer2.
Reads will also be discarded if they fail (minpolymer AND (ldf OR entropy)).
A read pair will be discarded if either read is discarded.

Depth-filtering parameters:
mincount=2      Minimum number of times a read kmer must occur in the 
                read set to be considered 'high-depth'.
ldf=0.24        (lowdepthfraction) Consider a read low-depth if at least
                this fraction of kmers are low depth.  Setting this above 1
                will disable depth analysis (making the program run faster).
ldf2=1.1        Discard reads with at least this fraction of low-depth kmers.
                Values above 1 disables this filter (e.g., for metagenomes).

Entropy-filtering parameters:
entropy=0.67    Reads with average entropy below this are considered 
                low-entropy.
entropy2=0.2    Reads with average entropy below this are discarded.

Quality-filtering parameters (only useful if q-scores are correct):
quality=12.5    Reads with average quality below this are considered 
                low-quality.
quality2=7.5    Reads with average quality below this are discarded.

Homopolymer-filtering parameters:
polymers=GC     Look for homopolymers of these symbols.  e.g., polymers=GC
                would look for poly-G or poly-C (but not poly-GC).
minpolymer=20   Minimum length of homopolymers.
minpolymer2=29  Discard any read with a homopolymer of at least this length.
purity=0.85     Min fraction of the homopolymer region that is the correct
                symbol.  For example, GGGGGGAGGG is length 10 with 9 Gs, for
                a purity of 0.90 (insufficient in this case due to length).

Trimming parameters:
trimpolymers=   Homopolymers to use for trimming.  If unspecified, it will
                be the same as 'polymers'.
trimleft=6      Trim left ends where there is a homopolymer at least this
                long; 0 disables trimming.
trimright=6     Trim left ends where there is a homopolymer at least this
                long; 0 disables trimming.
trim=           Sets both trimleft and trimright.
maxnonpoly=2    Trim through up to this many consecutive mismatches.
minlen=50       Discard reads shorter than this after trimming.

Other parameters:
quantize=1      If greater than 1, bin the quality scores to reduce file size.
cardinality=t   Report estimated number of unique output kmers.

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

polyfilter() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP bloom.PolyFilter $@"
	echo $CMD >&2
	eval $CMD
}

polyfilter "$@"
