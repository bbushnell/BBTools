#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 16, 2020

Description:  Provides some of Reformat's functionality in a ZMW-aware tool.

Usage:  reformatpb.sh in=<input file> out=<output file> outb=<bad reads>

File I/O parameters:
in=<file>       Primary input.
out=<file>      (outgood) Output for good reads.
outb=<file>     (outbad) Output for discarded reads.
stats=<file>    Print screen output here instead of to the screen.
json=f          Print stats as json.
schist=<file>   Subread count per ZMW histogram.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
kzt=f           (keepzmwstogether) Send all reads from a ZMW to the same file,
                either good or bad output.
minlen=40       Do not output reads shorter than this, after trimming.
ccsin=f         Input reads are CCS, meaning they are all full-pass.
                Not currently used for anything.
trimpolya=f     Trim terminal poly-A and poly-T sequences, for some isoseq
                libraries.
minpolymer=5    Don't trim poly-A sequence shorter than this.
polyerror=0.2   Max error rate for trimming poly-A.
flaglongreads=f    True to flag reads longer than 1.5x median to be discarded.
longreadmult=1.5   Multiplier to consider a read suspiciously long.

Whitelist and Blacklist Parameters:
whitelist=      ZMW identifiers, as a comma-delimited list of integers,
                or files with one integer per line.  All ZMWs not in the
                list will be discarded.
blacklist=      All ZMWs in the list will be discarded.

Sampling parameters (avoid using more than one of these at a time):
reads=-1        If positive, quit after processing this many reads.
zmws=-1         If positive, quit after processing this many ZMWs.
bestpass=f      Set to true to keep only the best read per ZMW.  This is
                the median length read of the non-outermost reads.
                If there are 2 or fewer passes, the longest will be chosen.
longestpass=f   Set to true to keep only the longest read per ZMW.
samplerate=1.0  Retain this fraction of input reads.
samplereadstarget=-1  If positive, retain this many reads.
samplebasestarget=-1  If positive, retain this many bases.
samplezmwstarget=-1   If positive, retain this many ZMWs.
subsamplefromends=f   If true, eliminate outermost reads first, then inner.

CCS Parameters (Note: CCS is still experimental)
ccs=f           Make a single consensus read per ZMW (CPU-intensive).
minpasses=0     Discard ZMWs with fewer than this many passes (estimated;
                first and last subread are usually partial).
minsubreads=0   Discard ZMWs with fewer than this many subreads.
reorient=f      Try aligning both strands in case ZMW ordering is broken.
minshredid=0.6  Do not include shreds with identity below this in consensus.

Entropy Parameters (recommended setting is 'entropy=t'):
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
#z3="-Xss16m"
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

reformatpb() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP icecream.ReformatPacBio $@"
	if [[ $silent != 1 ]]; then
		echo $CMD >&2
	fi
	eval $CMD
}

reformatpb "$@"
