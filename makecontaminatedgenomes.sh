#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 29, 2017

Description:  Generates synthetic contaminated partial genomes from clean genomes.
Output is formatted as (prefix)_bases1_fname1_bases2_fname2_counter_(suffix).

Usage:        makecontaminatedgenomes.sh in=<file> out=<pattern>

I/O parameters:
in=<file>       A file containing one input file path per line.
out=<pattern>   A file name containing a # symbol (or other regex).
                The regex will be replaced by source filenames.

Processing Parameters:
count=1         Number of output files to make.
seed=-1         RNG seed; negative for a random seed.
exp1=1          Exponent for genome 1 size fraction.
exp2=1          Exponent for genome 2 size fraction.
subrate=0       Rate to add substitutions to new genomes (0-1).
indelrate=0     Rate to add substitutions to new genomes (0-1).
regex=#         Use this substitution regex for replacement.
delimiter=_     Use this delimiter in the new file names.

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
}
calcXmx "$@"

makecontaminatedgenomes() {
	local CMD="java $EA $EOOM $z -cp $CP synth.MakeContaminatedGenomes $@"
	echo $CMD >&2
	eval $CMD
}

makecontaminatedgenomes "$@"
