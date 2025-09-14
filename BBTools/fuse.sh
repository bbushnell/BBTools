#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 11, 2023

Description:  Fuses sequences together, padding gaps with Ns.

Usage:   fuse.sh in=<input file> out=<output file> pad=<number of Ns>

Parameters:
in=<file>       The 'in=' flag is needed if the input file is not the 
                first parameter.  'in=stdin' will pipe from standard in.
out=<file>      The 'out=' flag is needed if the output file is not the 
                second parameter.  'out=stdout' will pipe to standard out.
pad=300         Pad this many N between sequences.
maxlen=2g       If positive, don't make fused sequences longer than this.
quality=30      Fake quality scores, if generating fastq from fasta.
overwrite=t     (ow) Set to false to force the program to abort rather 
                than overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change 
                compression level; lower compression is faster.
fusepairs=f     Default mode fuses all sequences into one long sequence.
                Setting fusepairs=t will instead fuse each pair together.
name=           Set name of output sequence.  Default is the name of
                the first input sequence.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding
                autodetection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

fuse() {
	local CMD="java $EA $EOOM $z -cp $CP synth.FuseSequence $@"
	echo $CMD >&2
	eval $CMD
}

fuse "$@"
