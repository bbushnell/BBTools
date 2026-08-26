#!/bin/bash

usage(){
echo "
Written by Neptune and Brian Bushnell
Last modified August 17, 2026

Description:  Greedy set-cover kmer selection for tRNA (or any short-gene)
covering sets.  Iteratively selects the most prevalent kmers from a pool
of sequences, evicts covered sequences, and repeats until a coverage
target or kmer budget is reached.

Two-tier ranking: each round selects the top 2*step kmers by current count
on the remaining pool, re-sorts by original prevalence (pre-eviction),
and keeps the top step.  This concentrates on kmers shared between rare
and common sequences rather than random sequence common only to stragglers.

Usage:	coveringset.sh in=<file> out=<file>

Input parameters:
in=<file>       Input FASTA of pool sequences (e.g. all tRNAs).
extra=<file>    Optional: additional sequences to weight the pool.
                Typically the consensus models; added copies= times.
copies=10       Number of copies of extra sequences to add to the pool.
rcomp=f         Canonicalize kmers: build the forward kmer and its reverse
                complement in the same rolling pass and count/select/match
                only the canonical (max) form, so a motif and its reverse
                complement are treated as one kmer instead of two.

Kmer parameters:
k=17            Kmer length for selection and output.
kdesign=        If set, select at this k but output at k.
                E.g. kdesign=18 k=17: select 18-mers, output their
                constituent 17-mers.  Each covered sequence then has
                >=2 kmer hits, enabling minkmerhits=2 for free.

Selection parameters:
step=500        Kmers to select per round.
maxkmers=       Maximum total kmers to select (0 = no limit).
target=0.999    Stop when this fraction of pool sequences is covered.

Java Parameters:
-Xmx            Set memory usage.  Default autodetected.
"
}

#This block allows symlinked shellscripts to correctly find the repo.
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

z="-Xmx1g"
z2="-Xms1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

coveringset() {
	local CMD="java $EA $EOOM $z $z2 --add-modules jdk.incubator.vector -cp $CP prok.CoveringSet $@"
	echo $CMD >&2
	eval $CMD
}

coveringset "$@"
