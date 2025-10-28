#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 6, 2015

Description:  Splits Nextera LMP libraries into subsets based on linker orientation:
LMP, fragment, unknown, and singleton.
Please read bbmap/docs/guides/SplitNexteraGuide.txt for more information.

Usage:  splitnextera.sh in=<file> out=<file> outf=<file> outu=<file> outs=<file>

For pairs in two files, use in1, in2, out1, out2, etc.

*** Note ***
For maximal speed, before running splitnextera, the linkers can be replaced with a constant first.

In other words, you can either do this (which is slightly faster):
bbduk.sh in=reads.fq out=replaced.fq ktmask=J k=19 hdist=1 mink=11 hdist2=0 literal=CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG
splitnextera.sh in=replaced.fq out=longmate.fq outf=frag.fq outu=unknown.fq outs=singleton.fq

Or this:
splitnextera.sh in=reads.fq out=longmate.fq outf=frag.fq outu=unknown.fq outs=singleton.fq mask=t


I/O Parameters:
in=<file>       Input reads.  Set to 'stdin.fq' to read from stdin.
out=<file>      Output for pairs with LMP orientation.
outf=<file>     Output for pairs with fragment orientation.
outu=<file>     Pairs with unknown orientation.
outs=<file>     Singleton output.
ow=f            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f           (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto       ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).

Processing Parameters:
mask=f          Set to true if you did not already convert junctions to some symbol, and it will be done automatically.
junction=J      Look for this symbol to designate the junction bases.
innerlmp=f      Generate long mate pairs from the inner pair also, when the junction is found in both reads.
rename=t        Rename read 2 of output when using single-ended input.
minlength=40    (ml) Do not output reads shorter than this.
merge=f         Attempt to merge overlapping reads before looking for junctions.
testmerge=0.0   If nonzero, only merge reads if at least the fraction of input reads are mergable.

Sampling Parameters:

reads=-1        Set to a positive number to only process this many INPUT reads (or pairs), then quit.
samplerate=1    Randomly output only this fraction of reads; 1 means sampling is disabled.
sampleseed=-1   Set to a positive number to use that prng seed for sampling (allowing deterministic sampling).

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

function splitnextera() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.SplitNexteraLMP $@"
	echo $CMD >&2
	eval $CMD
}

splitnextera "$@"
