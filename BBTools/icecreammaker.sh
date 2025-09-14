#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2020

Description:  Generates synthetic PacBio reads to mimic the chimeric
inverted repeats from 'triangle reads', aka 'ice cream cones' - 
reads missing one adapter.

Usage:  icecreammaker.sh in=<file> out=<file> reads=100k minlen=500 maxlen=5000

Standard parameters:
in=<file>       A reference genome fasta (optional).
out=<file>      Synthetic read output.
idhist=<file>   Identity histogram output.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Length parameters:
NOTE: "length" parameters dictate subread length (for normal reads).
"Movie" parameters dictate sequence length (with concatenated subreads).

minlen=500      (minlength) Minimum length of genomic sequence molecules.
maxlen=5000     (maxlength) Maximum length of genomic sequence molecules.
len=            (length) Set minlen and maxlen to the same number.
minmovie=500    (minmov) Minimum length of movies.
maxmovie=40k    (maxmov) Maximum length of movies.
movie=          (mov) Set minmov and maxmov to the same number.

Ice cream parameters:
missingrate=0   (missing) Fraction of reads missing an adapter.
hiddenrate=0    (hidden) Fraction of adapters not detected.
bothends=f      Allow missing or hiddden adapters on both ends.

Other parameters:
zmws            (reads) Number of ZMWs to generate.  There are actually
                multiple subreads per zmw.
ccs=f           Make CCS reads (one read per ZMW, full pass only).
                You still need to specify the error rate.
gc=0.6          If a random genome is generated, use this GC fraction.
genomesize=10m  If a random genome is generated, make it this big.
irrate=0.0      Add inverted repeats until this fraction of the genome
                is inverted repeats.
irminlen=500    Minimum length of inverted repeats.
irmaxlen=5000   Maximum length of inverted repeats
irlen=          Set minirlen and maxirlen to the same number.
miner=0.05      (minerrorrate) Minimum error rate.
maxer=0.28      (maxerrorrate) Maximum error rate.
er=             (errorrate) Set minerrorrate and maxerrorrate.
NOTE: You can alternatively set minid, maxid, or id.


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

z="-Xmx2g"
z2="-Xms2g"
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

icecreammaker() {
	local CMD="java $EA $EOOM $z -cp $CP icecream.IceCreamMaker $@"
	echo $CMD >&2
	eval $CMD
}

icecreammaker "$@"
