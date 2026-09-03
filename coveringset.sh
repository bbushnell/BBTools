#!/bin/bash

usage(){
echo "
Written by Neptune, Sayu, and Brian Bushnell
Last modified September 2, 2026

Description:  Greedy set-cover kmer selection for tRNA (or any short-gene)
covering sets.  Iteratively selects the most prevalent kmers from a pool
of sequences, evicts covered sequences, and repeats until a coverage
target or kmer budget is reached.  Also runs on amino-acid sequences in
the full or a reduced alphabet (alphabet=, key=), and on many protein
families at once (families=), one covering set per family.

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

Alphabet parameters (protein mode; default is nucleotide):
alphabet=nt     nt (default; 2-bit, rcomp allowed), amino (the 20 amino
                acids, 5-bit), or an explicit symbol list such as
                alphabet=ACDFGHIW (one class per letter).  Residues outside
                the alphabet (X, B, Z, *, gaps) reset the rolling kmer.
key=            Translation key applied on the fly (no re-encoded FASTA):
                slash-separated groups, e.g. key=AST/C/DEKNQR/FY/GP/H/ILMV/W
                (each group collapses to its first letter), or a preset
                name: legacy (amino8), c6, c7, c8, c9, c12, c14, mi7.
                Output kmers are written in the class letters.

Family-batch parameters (protein mode):
families=       A directory of per-family protein FASTAs (family id = file
                name without .faa/.fa/.fasta[.gz]) or a manifest TSV of
                family_id<TAB>path.  One covering set per family, all in one
                JVM (families in parallel).  Requires summary=.
summary=<file>  Per-family summary TSV: members, excluded, selected kmers,
                rounds, coverage, uncovered.  out= then receives the sets
                TSV (family_id, kmer, selection_rank, original_count, round)
                with #alphabet/#key/#k header lines.
exclude=<file>  Member ids (column 1) removed from every family pool before
                selection, e.g. a held-out query set.
minhits=1       Evict a sequence only once it contains this many DISTINCT
                selected kmers (minhits=2 guarantees two independent seeds
                per covered member; kdesign=k+1 is the cheaper alternative
                but its two kmers overlap).
maxfamilies=0   If >0, a kmer present in more than this many families' pools
                is ineligible for selection (cross-family specificity).
commit=         Free-text provenance written to the output headers.

Kmer parameters:
k=17            Kmer length for selection and output.
kdesign=        If set, select at this k but output at k.
                E.g. kdesign=18 k=17: select 18-mers, output their
                constituent 17-mers.  Each covered sequence then has
                >=2 kmer hits, enabling minkmerhits=2 for free.

Selection parameters:
step=500        Kmers to select per round.
stepfraction=0  Adaptive sizing for rounds 2+: target this fraction of the
                remaining sequences, projected from the previous round's
                recovered-sequences-per-kmer rate. 0 disables adaptation.
minstepmult=0.1 Lower bound on adaptive step, as a multiple of initial step.
maxstepmult=2   Upper bound on adaptive step, as a multiple of initial step.
step2boost=1    Multiplier applied only to the projected second-round step.
maxkmers=       Maximum total kmers to select (0 = no limit).
target=0.999    Stop when this fraction of pool sequences is covered.

Performance parameters (kmer counting is multithreaded via a partitioned
counter -- each thread buffers kmers per-partition and only locks a partition
briefly to flush a batch, so there's no serial merge step and no single-table
size ceiling; useful for large corpora, e.g. whole-domain reference sets):
partitions=     Number of counter shards. Default: nearest odd value >=
                max(15, 2*threads) -- odd avoids aliasing against the
                structured low bits of 2-bit-packed kmers. More shards means
                less lock contention between threads and a higher
                total-distinct-kmer ceiling (each shard is its own table),
                at the cost of a little fixed memory overhead.
bufsize=200     Kmers a thread buffers per partition before flushing (taking
                that partition's lock). Larger reduces lock frequency.

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
