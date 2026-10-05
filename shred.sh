#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 4, 2026
Description:  Shreds sequences into shorter, possibly overlapping sequences.

Usage: shred.sh in=<file> out=<file> length=<int>
       shred.sh in=<assembly.fa> out=<pieces.fa> k=<int>

File Parameters:
in=<file>       Input sequences.
out=<file>      Destination of output shreds.

Processing Parameters:
k=0            A positive k selects de Bruijn branch shredding instead of length
                shredding. Count canonical k-mers over the input, then cut
                after each k-mer with multiple present predecessors or successors.
                Repeat count alone does not cause a cut. Keeps every input base,
                including Ns and short tails; consecutive branches can give 1-bp
                pieces. Uses two passes over a regular FASTA/FASTQ input file.
                Incompatible with explicit length/minlen/maxlen/median/variance/
                mode flags, equal=t, nonzero overlap/increment, or maxns filtering.
                k=0 retains the ordinary length-based behavior below.
hashonly=t     For branch mode with k>31, store two hashes per k-mer instead of
                its full packed sequence. Hash keys use 16 bytes plus a 4-byte
                count per table slot, independent of k. Set false for full keys.
                Small-k tables and ordinary length shredding are unaffected.
length=500      Desired length of shreds if a uniform length is desired.
minlen=-1       Shortest allowed shred.  The last shred of each input sequence
                may be shorter than desired length if this is not set.
maxlen=-1       Longest shred length.  If minlength and maxlength are both
                set, shreds will use a random flat length distribution.
median=-1       Alternatively, setting median and variance will override
                minlen and maxlen.
variance=-1
linear          When maxlen is greater than minlen, the distribution can
                be linear, exp, or log (pick one as a flag).
overlap=0       Amount of overlap between successive shreds.
reads=-1        If nonnegative, stop after this many input sequences.
equal=f         Shred each sequence into subsequences of equal size of at most
                'length', instead of a fixed size.
qfake=30        Quality score, if using fastq output.
filetid=f       Name shreds with a tid parsed from the filename (e.g. tid_5).
headertid=f     Name shreds with a tid parsed from sequence headers.

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

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4000m" "--xms=4000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP synth.Shred $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" synth.Shred "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
