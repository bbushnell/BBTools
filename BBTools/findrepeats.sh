#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 5, 2023

Description:  Finds repeats in a genome.

No alignment is done; a sequence is considered to be a repeat of depth D if
all kmers within it have a depth of at least D.  Gaps of up to length G
consecutive kmers with lower counts may be allowed, which typically finds
far more and substantially longer repeats even with a small G.

Usage:  findrepeats.sh in=<input file> out=<output file>

Standard parameters:
in=<file>       Primary input (the genome fasta).
out=<file>      Primary output (list of repeats as TSV).  If no file is
                given this will be printed to screen; to suppress printing,
                use 'out=null'.
outs=<file>     (outsequence) Optional sequence output, for printing or 
                masking repeats.
overwrite=f     (ow) False ('f') forces the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
k=31            Kmer length to use (range is 1-31; 31 is recommended).
gap=0           Maximum allowed gap length within a repeat, in kmers.
                This allows inexact repeats.  Note that a 1bp mismatch
                will spawn a length K gap.
qhdist=0        Hamming distance within a kmer; allows inexact repeats.
                Values above 0 become exponentially slower.
minrepeat=0     Minimum repeat length to report, in bases.  Nothing shorter
                than kmer length will be found regardless of this value.
mindepth=2      Ignore copy counts below mindepth. 
maxdepth=-1     If positive, copy counts greater than this will be reported
                as this number.  This can greatly increase speed in rare 
                cases of thousand+ copy repeats.
preview=27      Print this many characters of the repeat per line.  Set to 0
                to suppress (may save memory).
mask=f          Write sequence with masked repeats to 'outs'. Possible values:
                   f: Do not mask.
                   t: Mask (by default, 't' or 'true' are the same as 'soft').
                   soft: Convert masked bases to lower case.
                   hard: Convert masked bases to 'N'.
                   Other characters: Convert masked bases to that character.
print=t         (printrepeats) Print repeat sequence to outs.  'print' and
                'mask' are mutually exclusive so enabling one will disable
                the other.
weak=f          (weaksubsumes) Ignore repeats that are weakly subsumed by 
                other repeats. A repeat is subsumed if there is another repeat
                with greater depth at the same coordinates.  Since a 3-copy 
                repeat is also a 2-copy repeat, only the 3-copy repeat will be
                reported.  However, in the case that the 3-copy repeat is
                inexact (has gaps) and the 2-copy repeat is perfect, both will
                be reported when 'weak=f' as per default.  If you set the
                'weak=t' flag, only the highest-depth version will be reported
                even if it has more gaps.  In either case all 3 repeats would
                be reported, but with 'weak=f' some copies would be reported
                twice for the same coordinates, once as a depth-2 perfect 
                repeat and again as a depth-3 imperfect repeat.

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

findrepeats() {
	local CMD="java $EA $EOOM $z -cp $CP repeat.RepeatFinder $@"
	echo $CMD >&2
	eval $CMD
}

findrepeats "$@"
