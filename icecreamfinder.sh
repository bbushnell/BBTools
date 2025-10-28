#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 6, 2020

Description:  Finds PacBio reads containing inverted repeats.
These are candidate triangle reads (ice cream cones).
Either ice cream cones only, or all inverted repeats, can be filtered.

Usage:  icecreamfinder.sh in=<input file> out=<output file> outb=<bad reads>

File I/O parameters:
in=<file>       Primary input.
out=<file>      (outgood) Output for good reads.
outa=<file>     (outambig) Output for with inverted repeats, but it is unclear
                whether that is natural or artifactual.
outb=<file>     (outbad) Output for reads suspected as chimeric.
outj=<file>     (outjunction) Output for junctions in inverted repeat reads.
stats=<file>    Print screen output here instead of to the screen.
json=f          Print stats as json.
asrhist=<file>  Adapter alignment score ratio histogram.
irsist=<file>   Inverted repeat alignment score ratio histogram.
ambig=          Determine where ambiguous reads are sent.  They will ALWAYS
                be sent to outa if specified.  If not, they will be sent
                to outg (good) unless overridden by this flag.  Options:
                   ambig=good:  Send ambiguous reads to outg.
                   ambig=bad:  Send ambiguous reads to outb.
                   ambig=good,bad:  Send ambiguous reads to outg and outb.
                   ambig=null:  Do not send to outg or outb.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
alignrc=t       Align the reverse-complement of the read to itself to look 
                for inverted repeats.
alignadapter=t  Align adapter sequence to reads.
adapter=        default: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
icecreamonly=t  (ico) Only remove suspected triangle reads.  Otherwise, all
                inverted repeats are removed.
ksr=t           (keepshortreads) Keep non-triangle reads from triangle ZMWs.
kzt=f           (keepzmwstogether) Send all reads from a ZMW to the same file.
targetqlen=352  (qlen) Make queries of this length from a read tip.
qlenfraction=0.15   Try to make queries at most this fraction of read length.
                For short reads this will override targetqlen.
minlen=40       Do not output reads shorter than this, after trimming.
minqlen=100     Do not make queries shorter than this.  For very short
                reads this will override qlenfraction.
shortfraction=0.4   Only declare a read to be a triangle if the short half 
                of the repeat is at least this fraction of read length.
ccs=f           Input reads are CCS, meaning they are all full-pass.
                In this case you should increase minratio.
trim=t          Trim adapter sequence from read tips.
trimpolya=f     Trim terminal poly-A and poly-T sequences, for some isoseq
                libraries.
minpolymer=5    Don't trim poly-A sequence shorter than this.
polyerror=0.2   Max error rate for trimming poly-A.


Speed and sensitivity parameters:
jni=f           Enable C code for higher speed and identical results.
minratio=       Fraction of maximal alignment score to consider as matching.
                Higher is more stringent; lower allows more sequencing errors.
                This is VERY SENSITIVE.  For error-corrected reads it should
                be set higher.  It is roughly the expected identity of one
                read to another (double the per-read error rate).
minratio1=0.59  Set minratio for the first alignment pass only.
minratio2=0.64  Set minratio for the second alignment pass only.
adapterratio=0.18   Initial adapter detection sensitivity; affects speed.
adapterratio2=.325  Final adapter detection sensitivity.
minscore=-800   Exit alignment early if score drops below this.

Entropy parameters (recommended setting is 'entropy=t'):
minentropy=-1   Set to 0.4 or above to remove low-entropy reads;
                range is 0-1, recommended value is 0.55.  0.7 is too high.
                Negative numbers disable this function.
entropyk=3      Kmer length for entropy calculation.
entropylen=350  Reads with entropy below cutoff for at least this many
                consecutive bases will be removed.
entropyfraction=0.5     Alternative minimum length for short reads; the shorter
                        of entropylen and entfraction*readlength will be used.
entropywindow=50        Window size used for entropy calculation.
maxmonomerfraction=0.74 (mmf) Also require this fraction of bases in each
                        window to be the same base.

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
JNI="-Djava.library.path=""$DIR""jni/"
#JNI=""

z="-Xmx2g"
z2="-Xms2g"
z3="-Xss16m"
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

icecream() {
	local CMD="java $EA $SIMD $EOOM $z $z2 $z3 $JNI -cp $CP icecream.IceCreamFinder $@"
	if [[ $silent != 1 ]]; then
		echo $CMD >&2
	fi
	eval $CMD
}

icecream "$@"
