#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 23, 2014

*** DEPRECATED: This should still work but is no longer maintained. ***

Description:  Annotates reads with their kmer depth.

Usage:        kmercoverage in=<input> out=<read output> hist=<histogram output>

Input parameters:
in2=null            Second input file for paired reads
extra=null          Additional files to use for input (generating hash table) but not for output
fastareadlen=2^31   Break up FASTA reads longer than this.  Can be useful when processing scaffolded genomes
tablereads=-1       Use at most this many reads when building the hashtable (-1 means all)
kmersample=1        Process every nth kmer, and skip the rest
readsample=1        Process every nth read, and skip the rest

Output parameters:
hist=null           Specify a file to output the depth histogram
histlen=10000       Max depth displayed on histogram
reads=-1            Only process this number of reads, then quit (-1 means all)
sampleoutput=t      Use sampling on output as well as input (not used if sample rates are 1)
printcoverage=f     Only print coverage information instead of reads
useheader=f         Append coverage info to the read's header
minmedian=0         Don't output reads with median coverage below this
minaverage=0        Don't output reads with average coverage below this
zerobin=f           Set to true if you want kmers with a count of 0 to go in the 0 bin instead of the 1 bin in histograms.
                    Default is false, to prevent confusion about how there can be 0-count kmers.
                    The reason is that based on the 'minq' and 'minprob' settings, some kmers may be excluded from the bloom filter.

Hashing parameters:
k=31                Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)
cbits=8             Bits per cell in bloom filter; must be 2, 4, 8, 16, or 32.  Maximum kmer depth recorded is 2^cbits.
                    Large values decrease accuracy for a fixed amount of memory.
hashes=4            Number of times a kmer is hashed.  Higher is slower.
                    Higher is MORE accurate if there is enough memory, and LESS accurate if there is not enough memory.
prefilter=f         True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable.
prehashes=2         Number of hashes for prefilter.
passes=1            More passes can sometimes increase accuracy by iteratively removing low-depth kmers
minq=7              Ignore kmers containing bases with quality below this
minprob=0.5         Ignore kmers with overall probability of correctness below this
threads=X           Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X by up to 4.

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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

kmercoverage() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.KmerCoverage prefilter=true bits=16 interleaved=false $@"
	echo $CMD >&2
	eval $CMD
}

kmercoverage "$@"
