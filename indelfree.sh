#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 6, 2026

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
outh=<file>     Sam header output (optional).  Due to the streaming nature,
                primary sam output is headerless, but this can be concatenated
		with the main sam file.
subs=5          (s) Maximum allowed substitutions.
minid=0.0       Minimum allowed identity.  Actual substitions allowed will be
                max(subs, (int)(qlen*(1-minid)))
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
iterations=200k Iterations for error distribution probability simulation.
qstep=1         Only look up every Nth query kmer (higher is faster).
rstep=1         Only index every Nth reference kmer.  Qstep is faster, but
                rstep uses less memory with long reference sequences.

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

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$0"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1000m" "--xms=1000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.IndelFreeAligner $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
