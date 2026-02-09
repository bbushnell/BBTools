#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 8, 2026

Description:  Aligns sequences, not allowing indels.
Brute force mode guarantees all alignments will be found and reported,
up to the maximum allowed number of substitutions.
Indexed mode uses an adaptive Multi-K strategy.  Queries are binned by
length and error rate, and the reference is indexed with multiple kmer
lengths (e.g. k=10,12,14) to optimize speed without sacrificing sensitivity.
This loads all reads into memory and streams the reference, unlike
a traditional aligner, so it is designed for a relatively small query set
and potentially enormous reference set.
Speed and sensitivity are greatly affected by the list of kmer lengths (k),
max subs allowed (subs), minimum identity (minid), minumum seed 
probability (minprob), and minimum query length.
Specifically, speed can be increased by:
Eliminating short values of k (such changing the default to k=10,12,14);
Decreasing minsubs;
Increasing minid;
Decreasing minprob;
Increasing minqlen.

Usage:  indelfree.sh in=spacers.fa ref=contigs.fa out=mapped.sam

Parameters:
in=<file>       Query input.  These will be stored in memory.
ref=<file>      Reference input.  These will be streamed.
out=<file>      Sam output.
outh=<file>     Sam header output (optional).  Due to the streaming nature,
                primary sam output is headerless, but this can be concatenated
                with the main sam file.
subs=5          (s) Maximum allowed substitutions.
minid=0.85      Minimum allowed identity.  Actual substitions allowed will be
                max(subs, (int)(qlen*(1-minid)))
minqlen=1       Ignore queries shorter than this.
minrlen=1       Ignore reference sequences shorter than this.
simd=t          Enable SIMD alignment.
threads=        Set the max number of threads; default is logical cores.
                Memory usage is proportional to threads times ref contig lengths.

Index Parameters:
index=t         If true, build a kmer index to accelerate search.
                Otherwise, brute force mode aligns queries to all locations.
k=8,9,10,12,14  Index kmer lengths (1-15).  Can be a single integer or a
                comma-delimited list.  The aligner will automatically select
                the longest valid K from the list for each query to maximize
                speed.  More lengths use more indexing time, but not more RAM.
		Short kmers (below 12) with short queries and high
		subs or low minid is very slow.
minprob=0.999   Calculate the number of seed hits needed, on a per-query
                basis, to ensure this probability of finding valid alignments.
                1 ensures optimality; 0 requires all seed hits; and negative
                numbers disable this, using the minhits setting only.
                When enabled, the min hits used for a query is the maximum
                of the minhits flag and the probabilistic model.
                This setting also controls the value of k chosen for a query.
mm=1            Middle mask length; the number of wildcard bases in the kmer.
                Must be shorter than k-1; 0 disables middle mask.
blacklist=2     Blacklist homopolymer kmers up to this repeat length.
chunk=1m        Fuse short sequences into chunks this long for indexing.
                Longer can be faster, but uses more memory.
minhits=1       Require this many seed hits to perform alignment.
prescan=t       Count query hits before filling seed location lists.
list=t          Store seed hits in lists rather than maps.
                Maps are optimized for shorter kmers and more positive hits.
iterations=200k Iterations for error distribution probability simulation.
qstep=1         Only look up every Nth query kmer (higher is faster).
rstep=1         Only index every Nth reference kmer.  Qstep is faster, but
                rstep uses less memory with long reference sequences.

Entropy parameters:
emask=f         Entropy-mask reference sequences to reduce low-complexity 
                spurious matches.
ewindow=80      Use this window length for entropy calculation.
ek=4            Kmer length for entropy calculation.
ecutoff=0.7     Mask windows with entropy below this.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ifa.IndelFreeAligner4 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
