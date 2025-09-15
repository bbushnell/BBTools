#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 22, 2023

Description:  Filters sequences by exact sequence matches.
This can also handle inexact matches but is extremely slow in that mode.

Usage:  filterbysequence.sh in=<file> out=<file> ref=<file> include=<t/f>

I/O Parameters:
in=             Primary input. 'in2' will specify a second file.
out=            Primary out. 'out2' will specify a second file.
ref=            A reference file or comma-delimited list of files.
literal=        A literal sequence or comma-delimited list of sequences.
ow=t            (overwrite) Overwrites files that already exist.

Processing Parameters:
include=f       Set to 'true' to include the filtered sequences rather
                than excluding them.
rcomp=t         Match reverse complements as well.
case=f          (casesensitive) Require matching case.
storebases=t    (sb) Store ref bases.  Requires more memory.  If false,
                case-sensitive matching cannot be done, and the matching
                will be probabilistic based on 128-bit hashcodes.
threads=auto    (t) Specify the number of worker threads.
subs=0          Maximum number of substitutions allowed.
mf=0.0          (mismatchfraction) Maximum fraction of bases that can
                mismatch.  The actual number allowed is the max of subs
                and mf*min(query.length, ref.length).
lengthdif=0     Maximum allowed length difference between query and ref.

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

z="-Xmx800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function filterbysequence() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.FilterBySequence $@"
	echo $CMD >&2
	eval $CMD
}

filterbysequence "$@"
