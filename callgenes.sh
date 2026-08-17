#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 16, 2026

Description:  Finds orfs and calls genes in unspliced prokaryotes.
This includes bacteria, archaea, viruses, and mitochondria.
Can also predict 16S, 18S, 23S, 5S, and tRNAs.

Usage:  callgenes.sh in=contigs.fa out=calls.gff outa=aminos.faa out16S=16S.fa

File parameters:
in=<file>       A fasta file; the only required parameter.
out=<file>      Output gff file.
outa=<file>     Amino acid output.
out16s=<file>   16S output.
model=<file>    A pgm file or comma-delimited list.
                If unspecified a default model will be used.
stats=stderr    Stats output (may be stderr, stdin, a file, or null).
hist=null       Gene length histogram.
compareto=      Optional reference gff file to compare with the gene calls.
                'auto' will name it based on the input file name.

Formatting parameters:
json=false      Print stats in JSON.
binlen=21       Histogram bin length.
bins=1000       Maximum histogram bins.
pz=f            (printzero) Print histogram lines with zero count.



Taxonomy parameters:
taxonomy=t      Use QuickClade to classify the input and select a phylum-
                specific gene model for improved accuracy.  Intended for
                ISOLATES (one organism per file): a metagenome would be
                misclassified to a single phylum, so set taxonomy=f for
                mixed samples.  Requires a QuickClade server; if it is
                unreachable, callgenes falls back to the general model with
                a warning (it never fails on this account).
percontig=f     Classify each contig separately (for metagenomes).
                Default is per-file (classify once for all contigs).
taxaddress=     QuickClade server address.  Default: refseq.

tRNA detection parameters:
scavengeonly=t  Call tRNAs with the kmer-guided scavenger only (finds tRNAs at
                conserved long-kmer positions).  This is the shipped default;
                set scavengeonly=f to use the PGM candidate generator instead.
scavenge=f      Run the scavenger IN ADDITION to the PGM candidate generator
                (augment mode) rather than replacing it.
trnaintron=t    Intron-aware pass: splice candidate spans before verification
                to recover intron-containing (mainly archaeal) tRNAs.
                (Not used while scavengeonly is on.)
mintrnakhits=1  Min conserved tRNA long-kmers a candidate must carry before it
                is aligned; a cheap pre-filter that rejects non-tRNA windows.
                0 disables the pre-filter.

tRNA alignment parameters:
trnaalign=t     Align predicted tRNAs to a consensus library to verify
                and annotate them.  Reduces false positives and adds
                anticodon/amino-acid annotation.  Uses built-in library
                by default; override with trnalib= and trnamodel=.
trnalib=<file>  Custom tRNA consensus library (fasta).
trnamodel=<file> Custom tRNA HBM model file.
indexk=7        Kmer length for the library shortlist index.  Longer is more
                selective (fewer models aligned per candidate); shorter is more
                permissive.  7 is the shipped default.
indextopn=60    Max library models aligned per candidate (search breadth).
indexminhits=12 Min shared index-kmers for a model to enter the shortlist.
                Fixed fallback; used only when adaptiveminhits=f.
adaptiveminhits=t  Adapt the shortlist cutoff per candidate instead of the
                fixed indexminhits (the shipped default).  Cutoff = ceil(max(
                adaptfloor, adapttopfrac*maxSharedKmers, adaptqfrac*queryKmers)).
adaptfloor=11   Absolute floor for the adaptive cutoff (the constant term).
adapttopfrac=0.48  Adaptive cutoff as a fraction of the best model's shared-
                kmer count.
adaptqfrac=0.072   Adaptive cutoff as a fraction of the candidate's kmer count.
patience=20     Stop aligning after this many models without improvement
                once a passing hit has been found (with earlyexit).
earlyexit=t     Enable the patience-based early exit.
idpass=0.75     Min alignment identity to accept a tRNA.
hbmpass=0.75    Min HBM model score to accept a borderline tRNA.
acextract=t     Extract the anticodon directly from each verified tRNA's
                structure (anticodon loop position projected through the
                alignment).  Adds an anticodon: attribute to the GFF.
acvalidate=12   Min structural score (anticodon stem + U33 + purine-37,
                max 15) to trust an extracted anticodon; failures fall
                back to model-name annotation.
acmargin=3      Min score margin between the best and runner-up anticodon
                register; ambiguous positions fall back to the model name.
trimtrna=t      Trim verified tRNA boundaries to the consensus alignment
                extent, removing scanner slop for accurate coordinates.
maxtrna=120     (Experimental) Raise the tRNA candidate length cap, enabling
                relaxed length scoring for over-length candidates such as
                intron-containing archaeal tRNAs.  Measured neutral: unspliced
                candidates still fail alignment verification.

Advanced tRNA candidate-generation thresholds (rarely changed):
trnaregion=20   Region-open score cutoff.
trnacand=36     Composite candidate score cutoff (tRNA is very sensitive here).
trnastart=2.4   Start point-model score cutoff.
trnastop=1.5    Stop point-model score cutoff.
trnainner=2.2   Average inner-kmer score cutoff.

Other parameters:
minlen=60       Don't call genes shorter than this.
trd=f           (trimreaddescription) Set to true to trim read headers after
                the first whitespace.  Necessary for IGV.
merge=f         For paired reads, merge before calling.
detranslate=f   Output canonical nucleotide sequences instead of amino acids.
recode=f        Re-encode nucleotide sequences over called genes, leaving
                non-coding regions unchanged.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

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

	parseJavaArgs "--xmx=6g" "--xms=6g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.CallGenes $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
