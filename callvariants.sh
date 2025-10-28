#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 31, 2025

Description:  Calls variants from sam or bam input.
In default mode, all input files are combined and treated as a single sample.
In multisample mode, each file is treated as an individual sample,
and gets its own column in the VCF file.  Unless overridden, input file
names are used as sample names.
Please read bbmap/docs/guides/CallVariantsGuide.txt for more information,
or bbmap/pipelines/variantPipeline.sh for a usage example.

Usage:  callvariants.sh in=<file,file,...> ref=<file> vcf=<file>

Input may be sorted or unsorted.
The reference should be fasta.

I/O parameters:
in=<file>       Input; may be one file or multiple comma-delimited files.
list=<file>     Optional text file containing one input file per line.
                Use list or in, but not both.
out=<file>      Output variant list in var format.  If the name ends
                with .vcf then it will be vcf format.
vcf=<file>      Output variant list in vcf format.
outgff=<file>   Output variant list in gff format.
ref=<file>      Reference fasta.  Required to display ref alleles.
                Variant calling will be more accurate with the reference.
vcfin=<file>    Force calls at these locations, even if allele count is 0.
shist=<file>    (scorehist) Output for variant score histogram.
zhist=<file>    (zygosityhist) Output for zygosity histogram.
qhist=<file>    (qualityhist) Output for variant base quality histogram.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
extended=t      Print additional variant statistics columns.
sample=         Optional comma-delimited list of sample names.
multisample=f   (multi) Set to true if there are multiple sam/bam files,
                and each should be tracked as an individual sample.
vcf0=           Optional comma-delimited list of per-sample outputs.
                Only used in multisample mode.
bgzip=t         Use bgzip for gzip compression.
samstreamer=t   (ss) Load reads multithreaded to increase speed.
                Disable to reduce the number of threads used.  The number of
                streamer threads can be set with e.g. 'ss=4'; default is 6.
streamermf=8    (ssmf) Allow multiple sam files to be read simultaneously.
                Set ssmf=X to specify the maximum number or ssmf=f
                to disable.

Processing Parameters:
prefilter=f     Use a Bloom filter to exclude variants seen fewer than
                minreads times.  Doubles the runtime but greatly reduces
                memory usage.  The results are identical.
coverage=t      (cc) Calculate coverage, to better call variants.
ploidy=1        Set the organism's ploidy.
rarity=1.0      Penalize the quality of variants with allele fraction 
                lower than this.  For example, if you are interested in
                4% frequency variants, you could set both rarity and
                minallelefraction to 0.04.  This is affected by ploidy - 
                a variant with frequency indicating at least one copy
                is never penalized.
covpenalty=0.8  (lowcoveragepenalty) A lower penalty will increase the 
                scores of low-coverage variants, and is useful for 
                low-coverage datasets.
useidentity=t   Include average read identity in score calculation.
usepairing=t    Include pairing rate in score calculation.
usebias=t       Include strand bias in score calculation.
useedist=t      Include read-end distance in score calculation.
homopolymer=t   Penalize scores of substitutions matching adjacent bases.
nscan=t         Consider the distance of a variant from contig ends when 
                calculating strand bias.
callsub=t       Call substitutions.
calldel=t       Call deletions.
callins=t       Call insertions.
calljunct=f     Call junctions (in development).
nopassdot=f     Use . as genotype for variations failing the filter.

Coverage Parameters (these mainly affect speed and memory use):
32bit=f         Set to true to allow coverage tracking over depth 65535,
                which increases memory use.  Variant calls are impacted
                where coverage exceeds the maximum.
atomic=auto     Increases multithreaded speed; forces 32bit to true.
                Defaults to true if there are more than 8 threads.
strandedcov=f   (stranded) Tracks per-strand ref coverage to print the MCOV
                and DP4 fields.  Requires more memory when enabled.  Strand
                of variant reads is tracked regardless of this flag.

Trimming Parameters:
border=5        Trim at least this many bases on both ends of reads.
qtrim=r         Quality-trim reads on this end
                   r: right, l: left, rl: both, f: don't quality-trim.
trimq=10        Quality-trim bases below this score.

Realignment Parameters:
realign=f       Realign all reads with more than a couple mismatches.
                Decreases speed.  Recommended for aligners other than BBMap.
unclip=f        Convert clip symbols from exceeding the ends of the
                realignment zone into matches and substitutitions.
repadding=70    Pad alignment by this much on each end.  Typically,
                longer is more accurate for long indels, but greatly
                reduces speed.
rerows=602      Use this many rows maximum for realignment.  Reads longer
                than this cannot be realigned.
recols=2000     Reads may not be aligned to reference seqments longer 
                than this.  Needs to be at least read length plus
                max deletion length plus twice padding.
msa=            Select the aligner.  Options:
                   MultiStateAligner11ts:     Default.
                   MultiStateAligner9PacBio:  Use for PacBio reads, or for
                   Illumina reads mapped to PacBio/Nanopore reads.

Sam-filtering Parameters:
minpos=         Ignore alignments not overlapping this range.
maxpos=         Ignore alignments not overlapping this range.
minreadmapq=4   Ignore alignments with lower mapq.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
secondary=f     Include secondary alignments.
supplementary=f Include supplementary alignments.
duplicate=f     Include reads flagged as duplicates.
invert=f        Invert sam filters.

Variant-Calling Cutoff Parameters:
minreads=2              (minad) Ignore variants seen in fewer reads.
maxreads=BIG            (maxad) Ignore variants seen in more reads.
mincov=0                Ignore variants in lower-coverage locations.
maxcov=BIG              Ignore variants in higher-coverage locations.
minqualitymax=15        Ignore variants with lower max base quality.
minedistmax=20          Ignore variants with lower max distance from read ends.
minmapqmax=0            Ignore variants with lower max mapq.
minidmax=0              Ignore variants with lower max read identity.
minpairingrate=0.1      Ignore variants with lower pairing rate.
minstrandratio=0.1      Ignore variants with lower plus/minus strand ratio.
minquality=12.0         Ignore variants with lower average base quality.
minedist=10.0           Ignore variants with lower average distance from ends.
minavgmapq=0.0          Ignore variants with lower average mapq.
minallelefraction=0.1   Ignore variants with lower allele fraction.  This
                        should be adjusted for high ploidies.
minid=0                 Ignore variants with lower average read identity.
minscore=20.0           Ignore variants with lower Phred-scaled score.
clearfilters            Clear all filters.  Filter flags placed after
                        the clearfilters flag will still be applied.

There are additionally max filters for score, quality, mapq, allelefraction,
and identity.

Other Parameters:
minvarcopies=0          If set to 0, a genotype (vcf GT field) of 0 or 0/0
                        will be called if observed allele frequency suggests
                        this is a minor allele.  If set to 1, GT field will
                        contain at least one 1.

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

callvariants() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP var2.CallVariants $@"
	echo $CMD >&2
	eval $CMD
}

callvariants "$@"
