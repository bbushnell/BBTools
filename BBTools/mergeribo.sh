#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 15, 2025

Description:  Merges files of SSU sequences to keep one per taxID.
By default, a consensus is generated per TaxID, then the sequence
best matching that consensus is used:
First, all sequences per TaxID are aligned to a reference consensus.
Second, the best-matching sequence is used as a seed, and all other
sequences for that TaxID are aligned to the seed to generate a new consensus.
Third, in 'consensus' mode, that consensus is simply output.
In 'best' mode (default), all sequences are aligned again to the new consensus,
and the best-matching is output.

Usage:  mergeribo.sh in=<file,file> out=<file>

Standard parameters:
in=<file,file>  Comma-delimited list of files.
out=<file>      Output file.
out2=<file>     Read 2 output if reads are in two files.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.
fastawrap=70    4000 is recommended to minimize filesize.

Processing parameters:
alt=<file>      Lower priority data.  Only used if there is no SSU associated
                with the TaxID from the primary input.
best=t          Output the best representative per taxID.
consensus=f     Output a consensus per taxID instead of the best input
                sequence.  Mutually exclusive with best.
fast=f          Output the best sequence based on alignment to global consensus
                (the seed) rather than individual consensus.
minid=0.62      Ignore sequences with identity lower than this to the global
                consensus.
maxns=-1        Ignore sequences with more than this many Ns, if non-negative.
minlen=1        Ignore sequences shorter than this.
maxlen=4000     Ignore sequences longer than this.
16S=t           Align to 16S consensus to pick the seed. Mutually exclusive.
18S=f           Align to 18S consensus to pick the seed. Mutually exclusive.
level=          If specified with a term like 'species' or 'genus', nodes
                will be promoted to that level, minimum, before consensus.
dada2=f         Output headers in dada2 format.

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
	freeRam 4000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

mergeribo() {
	local CMD="java $EA $EOOM $z -cp $CP prok.MergeRibo $@"
	echo $CMD >&2
	eval $CMD
}

mergeribo "$@"
