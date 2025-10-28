#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 26, 2021

Description:  Given a reference, generates a kmer spectrum including a 
specified number of substitutions, insertions, and deletions.  The output is 
useful for analyzing barcodes or other short oligos, and filtering using
BBDuk or Seal. Input may be fasta or fastq, compressed or raw.
See also kcompress, kmercountexact, and bbkmerset.

Usage:  kmutate.sh in=<file> out=<file> k=<kmer length> edist=<edit distance>

Examples:

        kmutate.sh in=x.fa out=y.fa k=31 edist=2
This will generate all 31-mers in x.fa, along with all 31-mer mutants with
an edit distance of up to 2.  For example, 1 substitution, or 1 substitution
plus 1 deletion, or any other combination of subsitutions, insertions,
and deletions that sum to 0, 1, or 2.

        kmutate.sh in=x.fa out=y.fa k=31 idist=1 ddist=3
This will generate all 31-mers in x.fa, along with all 31-mer mutants allowing
up to 1 insertion and 3 deletions, but no substitutions.  For example,
1 insertion and 3 deletions is possible (edit distance 4), but 1 deletion and 
1 substitution is not directly possible (though some equivalent mutants would
still be generated because a deletion and subsequent insertion is equivalent 
to a substitution).

Note that deletion events have limitations; e.g., they cannot occur on input
sequences of length k because the resulting sequence is shorter than k.  As
a result, running the program twice consecutively with edist=1 will not give
the same results as running once with edist=2 since the intermediate file will
be only k-length sequences.  However, running consecutively with sdist or
idist would be equivalent since they do not involve deletions.

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary output, or read 1 output.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
k=31            Kmer length; 1-31.
rcomp=t         Consider kmers equivalent to their reverse-complements.

Edit mode parameters (used if edist>0):
edist=0         Set the maximal edit distance (0-3).
smax=99         (optional) Don't allow more than this many total substitutions.
dmax=99         (optional) Don't allow more than this many total deletions.
imax=99         (optional) Don't allow more than this many total insertions.

SDI mode parameters:
sdist=0         Maximum substitutions allowed.
idist=0         Maximum insertions allowed.
ddist=0         Maximum deletions allowed (0-3).
emax=99         (optional) Don't allow more than this many total edits.

*** Note - please use SDI mode flags OR Edit mode flag, not both. ***
*** Both modes are equivalent, they just have different defaults. ***

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

kmutate() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP jgi.KExpand $@"
	echo $CMD >&2
	eval $CMD
}

kmutate "$@"
