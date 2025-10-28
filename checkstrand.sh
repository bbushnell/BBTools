#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 11, 2025

Description:  Estimates the strandedness of a library without alignment; 
intended for RNA-seq data.  Only the reads are required input to determine
strandedness, so this works when run with just a fastq file.  If sam/bam
input is used, additional alignment-based metrics will be reported.
If an assembly and gff file are provided, the affinity of reads to the plus 
or minus (sense or antisense) strand will also be calculated.  If a genome
is specified with no gff file, the genes will be called automatically with a
prokaryotic gene caller.  Strandedness and P/(P+M) ratios are similar but 
calculated in different ways, and thus will not exactly agree, but should be 
close.  For most calculations, only read 1 is used, or the merge of read 1
and read 2 if the merge flag is enabled and they overlap.

Usage:  checkstrand.sh in=<input file>


Output meaning:

Depth_Analysis: Based on comparing the fraction of kmers in forward and
                reverse orientation to a binomial distribution.
Strandedness:   Percent of reads that came from the majority strand, based
                on kmer depth.
StrandednessN:  Depth-normalized strandedness, where each unique kmer
                contributes equally regardless of depth.
AvgKmerDepth:   Average depth of kmers; typically higher than read depth.
Kmers>Depth1:   Fraction of kmers with depth over 1.  Singleton kmers cannot
                be used to calculate strandedness from depth.

Stop_Codon_Analysis:  Based on counting stop codons.
MajorStrandORF: Predicted major strand based on stop-codon analysis.
AvgReadLen:     Average length of R1, or the merged read pair.
AvgORFLen+:     Average ORF length in the best frame on the plus side of R1.
AvgORFLen-:     Average ORF length in the best frame on the minus side of R1.
AvgStopCount+:  Average number of stop codons in the best frame on plus side.
AvgStopCount-:  Average number of stop codons in the best frame on minus side.
GC_Content:     GC content of reads, which affects stop codon frequency.

PolyA_Analysis: Based on counting poly-A or poly-T tails.
MajorStrandPA:  Predicted major strand based on poly-A tail analysis.
                This is unreliable if the poly-A fraction is low.
PolyA/(PA+PT):  Ratio of (reads with poly-A tails)/(poly-A + poly-T tails).
PolyAFraction:  Fraction of reads ending in poly-A or poly-T.

Ref_Analysis:   Compares read kmer frequencies to a reference (if present).
                The reference can be a transcriptome or genome, and a gff file
                will be used if provided.
MajorStrandREF: The strand containing a majority of forward read kmers.
P/(P+M)_Ratio:  P is the sum of counts of plus kmers, and M is minus kmers.
GeneCoverage:   Fraction of transcriptome kmers represented in reads.
GenePrecision:  Fraction of read kmers found in the transcriptome.

Read_Gene_Calling_Analysis:  Uses gene-calling on reads to calculate which
                strand better fits a prokaryotic gene model.
MajorStrandRGC: Predicted major strand based on read gene calling.
P/(P+M)_Ratio:  P is the read count best matching the plus strand; M is minus.
AvgScorePlus:   Average score of called plus-strand genes.
AvgScoreMinus:  Average score of called plus-strand genes.
UsedFraction:   Fraction of reads with any called genes (or partial genes);
                this can be increased by merging the reads for longer frames.

Alignment_Results:  Requires sam/bam input.  The reads must have been
                mapped to a transcriptome or RNA-seq assembly, or to a 
                specified genome, or a gff file must be provided.
StrandednessAL: Percent of reads aligned to the dominant strand.  More 
                accurate for transcriptome-mapped than genome-mapped reads.
StrandednessAN: Depth-normalized strandedness, where each feature or
                contig contributes equally.
MajorStrandAL:  Strand to which a majority of reads aligned.
P/(P+M)_Ratio:  P is the number of plus-mapped reads, M is minus.
P/(P+M)_RatioN: Depth-normalized plus/total ratio.
PlusFeatures:   Fraction of features with majority plus-mapped reads.
AlignmentRate:  Fraction of reads that aligned.
Feature-Mapped: Fraction of reads that aligned to a feature in the gff.


Running on a fastq is simple, but there are multiple ways to run CheckStrand
on aligned data (in=, ref=, and gff= flags are not needed if the files have
proper extensions):

#1) This won't give alignment results, just kmer results
checkstrand.sh mapped.sam

#2) This will do gene-calling and the alignment strandedness will be based 
    on gene sense strand, but only works for prokaryotes/viruses
checkstrand.sh mapped.sam contigs.fa

#3) This will use the annotation and the alignment strandedness will be based
    on gene sense strand, works for proks, and should work for eukaryotes 
    (there are lots of ways to annotate multi-exon genes though)
checkstrand.sh mapped.sam genes.gff

#4) This will assume that the reference was a sense-strand transcriptome,
    and the alignment strandedness will be based on contig plus strand
checkstrand.sh mapped.sam transcriptome

#5) This will assume that the reference was unstranded contigs assembled
    from RNA-seq data, so the alignment strandedness will be based on
    contig majority strand
checkstrand.sh mapped.sam rnacontigs


Standard parameters:
in=<file>       Primary input (a fastq, fasta, sam or bam file).
in2=<file>      Secondary input for paired fastq in twin files.  Read 2 is 
                ignored unless merge=t.
out=<file>      Optional destination to redirect results instead of stdout.
outp=<file>     Optional output for plus-mapped reads.     
outm=<file>     Optional output for minus-mapped reads.
*Note: outp/outm require sam/bam input and either transcriptome mode or a gff.
The destination of plus-mapped r1 would be outp, but outm for plus-mapped r2.

Processing parameters:
ref=<file>      Optional reference (assembly) input.
gff=<file>      Optional gene annotation file input.
scafreport=<file>  Optional per-scaffold strandedness output.
transcriptome=f Set this to 't' if the reference is a sense-strand 
                transcriptome (rather than a genome assembly).  This applies
                to either a reference specified by 'ref' or the reference
                used for alignment, fo sam/bam input.
rnacontigs=f    Set this to 't' if the reference is contigs assembled from
                RNA-seq data, but with unknown orientation.  Only affects
                alignment results.
size=80000      Sketch size; larger may be more precise.
merge=f         Attempt to merge paired reads, and use the merged read when
                successful.  If unsuccessful only R1 is used.  This has a 
                performance impact on CPUs with few cores.
orf=t           Analyze stop codons and open reading frames.  Usually this
                will allow major strand determination without a reference, 
                unless GC is very high.
callreads=t     Perform gene-calling on reads using a prokarotic gene model.
                Not very relevant to eukaryotes.
passes=2        Two passes refines the gene model for better gene-calling on
                reads.  Only used if there is a reference.
samplerate=1.0  Set to a lower number to subsample the input reads; increases
                speed on CPUs with few cores.
sampleseed=17   Positive numbers are deterministic; negative use random seeds.
reads=-1        If positive, quit after processing this many reads.

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

z="-Xmx2g"
#z2="-Xms2g"
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

checkstrand() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.CheckStrand2 $@"
	echo $CMD >&2
	eval $CMD
}

checkstrand "$@"
