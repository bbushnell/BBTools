#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 15, 2025

Description:  Aligns sequences, not allowing indels.
Brute force mode guarantees all alignments will be found and reported,
up to the maximum allowed number of substitutions.
Indexed mode may remove this guarantee (depending on kmer length,
query length, and number of substitutions) but can be much faster.
This loads all reads into memory and streams the reference, unlike
a traditional aligner, so it is designed for a relatively small query set
and potentially enormous reference set.

Usage:  indelfree.sh in=spacers.fa ref=contigs.fa out=mapped.sam

Parameters:
in=<file>       Query input.  These will be stored in memory.
ref=<file>      Reference input.  These will be streamed.
out=<file>      Sam output.
subs=5          Maximum allowed substitutions.
simd            Enable SIMD alignment.  Only accelerates brute force mode.
threads=        Set the max number of threads; default is logical cores.

Index Parameters:
index=t         If true, build a kmer index to accelerate search.
k=13            Index kmer length (1-15); longer is faster but less sensitive.
                Very short kmers are slower than brute force mode.
mm=1            Middle mask length; the number of wildcard bases in the kmer.
                Must be shorter than k-1; 0 disables middle mask.
blacklist=2     Blacklist homopolymer kmers up to this repeat length.
step=1          Only use every Nth query kmer.
minhits=1       Require this many seed hits to perform alignment.
minprob=0.9999  Calculate the number of seed hits needed, on a per-query
                basis, to ensure this probability of finding valid alignments.
                1 ensures optimality; 0 requires all seed hits; and negative
                numbers disable this, using the minhits setting only.
                When enabled, the min hits used for a query is the maximum
                of minhits and the probabilistic model.
prescan=t       Count query hits before filling seed location lists.
list=t          Store seed hits in lists rather than maps.
                Maps are optimized for shorter kmers and more positive hits.


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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    # Parse Java arguments with tool-specific defaults
    # Use auto mode with 84% of available RAM, minimum 1000MB
    parseJavaArgs "--mem=1000m" "--percent=84" "--mode=auto" "$@"
    
    # Set environment paths
    setEnvironment
    
    # Set the Java memory parameters
    z="-Xmx${RAM}m"
    z2="-Xms${RAM}m"
}
calcXmx "$@"

align() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.IndelFreeAligner $@"
	echo $CMD >&2
	eval $CMD
}

align "$@"
