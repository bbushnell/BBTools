#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 25, 2025

Description:  Creates a mutant version of a genome.
Also produces a VCF listing the added mutations.

Usage:  mutate.sh in=<input file> out=<output file> id=<identity>

I/O parameters:
in=<file>       Input genome.
out=<file>      Output mutant genome.
vcf=<file>      Output VCF file showing variations added.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
subrate=0       Substitution rate, 0 to 1.
insrate=0       Insertion rate, 0 to 1.
delrate=0       Deletion rate, 0 to 1.
indelrate=0     Sets ins and del rate each to half of this value.
maxindel=1      Max indel length.
indelspacing=3  Minimum distance between subsequent indels.
id=1            Target identity, 0 to 1; 1 means 100%.
                If this is used it will override subrate and indelrate;
                99% of the mutations will be substitutions, and 1% indels.
fraction=1      Genome fraction, 0 to 1; 1 means 100%.  A lower fraction
                will randomly select that fraction on a per-sequence basis,
                possibly incurring one chimeric junction per sequence.
                Not compatible with VCF output.
period=-1       If positive, place exactly one mutation every X bases.
prefix=         Set this flag to rename the new contigs with this prefix
                and a number.
amino=f         Treat the input as amino acid sequence.
ploidy=1        Set the ploidy.  ploidy>1 allows heterozygous mutations.
                This will create one copy of each input sequence per ploidy.
hetrate=0.5     If polyploid, fraction of mutations that are heterozygous.
nohomopolymers=f  If true, prevent indels in homopolymers that lead to
                ambiguous variant calls.  For example, inserting A between
                AC or deleting T from TTTT.  This is mainly for grading 
                purposes.  It does not fully solve the problem, but greatly
                improves concordance (reducing disagreements by 70%).
                NOTE! nohomopolymers is temporarily disabled.
pad=0           Add this many random bases to the ends of input sequences.
                Padleft and padright may also be specified independently.

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

mutate() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP synth.MutateGenome $@"
	echo $CMD >&2
	eval $CMD
}

mutate "$@"
