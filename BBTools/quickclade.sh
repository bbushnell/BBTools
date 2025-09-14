#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Assigns taxonomy to query sequences by comparing kmer
frequencies to those in a reference database.  Developed for taxonomic
assignment of metagenomic bins, but it can also run on a per-sequence basis.
QuickClade is extremely fast and uses little memory.  However, the accuracy
declines for incomplete genomes.  The recommended minimum sequence length
is not yet known, but lower values of k5dif are more likely to be correct
to a lower taxonomic level.  k5dif represents the sum of the absolute values
of the differences between the 5-mer frequency spectra, so the range is 0-1.
Because no marker genes are used, QuickClade should perform similarly for any
clade in the reference dataset.
While the default reference is taxonomically labeled, you can use whatever
you want as a reference, with or without taxonomic labels.

Usage Examples:
quickclade.sh query1.fa query2.fa query3.fa
or
quickclade.sh bins
or
quickclade.sh contigs.fa percontig out=results.tsv usetree


File Parameters:
in=<file,file>  Query files or directories.  Loose file or directory names are
                also permitted.  Input can be fasta, fastq, or spectra files;
                spectra files are made by cladeloader.sh.
ref=<file,file> Reference files; the current default is:
                /clusterfs/jgi/groups/gentech/homes/bbushnell/clade/refseq_main.spectra.gz
                It is plaintext, human-readable, and pretty small.
out=stdout      Set to a file to redirect output.  Only the query results will
                be written here; progress messages will still go to stderr.

Basic Parameters:
percontig       Run one query per contig instead of per file.
minlen=0        Ignore sequences shorter than this in percontig mode.
hits=1          Print this many top hits per query.
steps=7         Only search up to this many GC intervals (of 0.01) away from
                the query GC.
oneline         Print results one line per query, tab-delimited.
callssu=f       Call 16S and 18S for alignment to reference SSU.
                This will affect the top hit ordering only if hits>1.

Advanced Parameters (mainly for benchmarking):
printmetrics    Output accuracy statistics; mainly useful for labeled data.
                Labeled data should have 'tid_1234' or similar in the header.
                Works best with 'usetree'.
printqtid       Print query TaxID.
banself         Ignore records with the same TaxID as the query.  Makes the
                program behave like that organism is not in the reference.
simd            Use vector instructions to accelerate comparisons.
maxk=5          Can be set to 4 or 3 to restrict kmer frequency comparisons
                to smaller kmers.  This may improve accuracy for small
                sequences/bins, but slightly reduces accuracy for large
                sequences/bins.
ccm=1.0         Threshold for using pentamers; lower is faster.
ccm2=1.5        Threshold for using tetramers.
gcdif=0.07      Initial maximum GC difference.
strdif=0.10     Initial maximum strandedness difference.
gcmult=0.5      Max GC difference as a fraction of best 5-mer difference.
strmult=1.2     Max strandedness difference as a fraction of best 5-mer diff.
ee=t            Early exit; increases speed.
entropy         Calculate entropy for queries.  Slow; negligible utility.
heap=1          Number of intermediate comparisons to store.
usetree         Load a taxonomic tree for better grading for labeled data.
aligner=quantum Options include ssa2, glocal, drifting, banded, crosscut.
Distance Metrics:
abs             Use absolute difference of kmer frequencies.
cos             Use 1-cosine similarity of kmer frequencies.
euc             Use Euclidian distance.
hel             Use Hellinger distance.
abscomp         GC-compensated version of abs (default).
Note:  The distance metric strongly impacts ccm, gcmult, and strmult.
       Defaults are optimized for abscomp.

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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    # Parse Java arguments with quickclade-specific defaults
    # For quickclade: fixed 4g memory, never autodetect
    parseJavaArgs "--mem=4g" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

quickclade() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP bin.CladeSearcher $@"
	echo $CMD >&2
	eval $CMD
}

quickclade "$@"
