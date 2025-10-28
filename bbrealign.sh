#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 8, 2018

Description:  Realigns mapped reads to a reference.

Usage:  bbrealign.sh in=<file> ref=<file> out=<file>

Input may be a sorted or unsorted sam or bam file.
The reference should be fasta.

I/O parameters:
in=<file>       Input reads.
out=<file>      Output reads.
ref=<file>      Reference fasta.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Trimming parameters:
border=0        Trim at least this many bases on both ends of reads.
qtrim=r         Quality-trim reads on this end
                   r: right, l: left, rl: both, f: don't quality-trim.
trimq=10        Quality-trim bases below this score.

Realignment parameters:
unclip=f        Convert clip symbols from exceeding the ends of the
                realignment zone into matches and substitutitions.
repadding=70    Pad alignment by this much on each end.  Typically,
                longer is more accurate for long indels, but greatly
                reduces speed.
rerows=602      Use this many rows maximum for realignment.  Reads longer
                than this cannot be realigned.
recols=2000     Reads may not be aligned to reference seqments longer 
                than this.  Needs to be at least read length plus
                max deletion length plus twice padding.
msa=            Select the aligner.  Options:
                   MultiStateAligner11ts:     Default.
                   MultiStateAligner9PacBio:  Use for PacBio reads, or for
                   Illumina reads mapped to PacBio/Nanopore reads.

Sam-filtering parameters:
minpos=         Ignore alignments not overlapping this range.
maxpos=         Ignore alignments not overlapping this range.
minreadmapq=4   Ignore alignments with lower mapq.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
secondary=f     Include secondary alignments.
supplementary=f Include supplementary alignments.
invert=f        Invert sam filters.

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

bbrealign() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP var2.Realign $@"
	echo $CMD >&2
	eval $CMD
}

bbrealign "$@"
