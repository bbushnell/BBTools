#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 19, 2017

Description:  Normalizes read depth based on kmer counts.
Can also error-correct, bin reads by kmer depth, and generate a kmer depth histogram.
However, Tadpole has superior error-correction to BBNorm.
Please read bbmap/docs/guides/BBNormGuide.txt for more information.

Usage:     bbnorm.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>

Input parameters:
in=null             Primary input.  Use in2 for paired reads in a second file
in2=null            Second input file for paired reads in two files
extra=null          Additional files to use for input (generating hash table) but not for output
fastareadlen=2^31   Break up FASTA reads longer than this.  Can be useful when processing scaffolded genomes
tablereads=-1       Use at most this many reads when building the hashtable (-1 means all)
kmersample=1        Process every nth kmer, and skip the rest
readsample=1        Process every nth read, and skip the rest
interleaved=auto    May be set to true or false to force the input read file to ovverride autodetection of the input file as paired interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Output parameters:
out=<file>          File for normalized or corrected reads.  Use out2 for paired reads in a second file
outt=<file>         (outtoss) File for reads that were excluded from primary output
reads=-1            Only process this number of reads, then quit (-1 means all)
sampleoutput=t      Use sampling on output as well as input (not used if sample rates are 1)
keepall=f           Set to true to keep all reads (e.g. if you just want error correction).
zerobin=f           Set to true if you want kmers with a count of 0 to go in the 0 bin instead of the 1 bin in histograms.
                    Default is false, to prevent confusion about how there can be 0-count kmers.
                    The reason is that based on the 'minq' and 'minprob' settings, some kmers may be excluded from the bloom filter.
tmpdir=$TMPDIR      This will specify a directory for temp files (only needed for multipass runs).  If null, they will be written to the output directory.
usetempdir=t        Allows enabling/disabling of temporary directory; if disabled, temp files will be written to the output directory.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
rename=f            Rename reads based on their kmer depth.

Hashing parameters:
k=31                Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)
bits=32             Bits per cell in bloom filter; must be 2, 4, 8, 16, or 32.  Maximum kmer depth recorded is 2^cbits.  Automatically reduced to 16 in 2-pass.
                    Large values decrease accuracy for a fixed amount of memory, so use the lowest number you can that will still capture highest-depth kmers.
hashes=3            Number of times each kmer is hashed and stored.  Higher is slower.
                    Higher is MORE accurate if there is enough memory, and LESS accurate if there is not enough memory.
prefilter=f         True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable.  The prefilter is more memory-efficient because it uses 2-bit cells.
prehashes=2         Number of hashes for prefilter.
prefilterbits=2     (pbits) Bits per cell in prefilter.
prefiltersize=0.35  Fraction of memory to allocate to prefilter.
buildpasses=1       More passes can sometimes increase accuracy by iteratively removing low-depth kmers
minq=6              Ignore kmers containing bases with quality below this
minprob=0.5         Ignore kmers with overall probability of correctness below this
threads=auto        (t) Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X due to I/O threads.
rdk=t               (removeduplicatekmers) When true, a kmer's count will only be incremented once per read pair, even if that kmer occurs more than once.

Normalization parameters:
fixspikes=f         (fs) Do a slower, high-precision bloom filter lookup of kmers that appear to have an abnormally high depth due to collisions.
target=100          (tgt) Target normalization depth.  NOTE:  All depth parameters control kmer depth, not read depth.
                    For kmer depth Dk, read depth Dr, read length R, and kmer size K:  Dr=Dk*(R/(R-K+1))
maxdepth=-1         (max) Reads will not be downsampled when below this depth, even if they are above the target depth.            
mindepth=5          (min) Kmers with depth below this number will not be included when calculating the depth of a read.
minkmers=15         (mgkpr) Reads must have at least this many kmers over min depth to be retained.  Aka 'mingoodkmersperread'.
percentile=54.0     (dp) Read depth is by default inferred from the 54th percentile of kmer depth, but this may be changed to any number 1-100.
uselowerdepth=t     (uld) For pairs, use the depth of the lower read as the depth proxy.
deterministic=t     (dr) Generate random numbers deterministically to ensure identical output between multiple runs.  May decrease speed with a huge number of threads.
passes=2            (p) 1 pass is the basic mode.  2 passes (default) allows greater accuracy, error detection, better contol of output depth.

Error detection parameters:
hdp=90.0            (highdepthpercentile) Position in sorted kmer depth array used as proxy of a read's high kmer depth.
ldp=25.0            (lowdepthpercentile) Position in sorted kmer depth array used as proxy of a read's low kmer depth.
tossbadreads=f      (tbr) Throw away reads detected as containing errors.
requirebothbad=f    (rbb) Only toss bad pairs if both reads are bad.
errordetectratio=125   (edr) Reads with a ratio of at least this much between their high and low depth kmers will be classified as error reads.
highthresh=12       (ht) Threshold for high kmer.  A high kmer at this or above are considered non-error.
lowthresh=3         (lt) Threshold for low kmer.  Kmers at this and below are always considered errors.

Error correction parameters:
ecc=f               Set to true to correct errors.  NOTE: Tadpole is now preferred for ecc as it does a better job.
ecclimit=3          Correct up to this many errors per read.  If more are detected, the read will remain unchanged.
errorcorrectratio=140  (ecr) Adjacent kmers with a depth ratio of at least this much between will be classified as an error.
echighthresh=22     (echt) Threshold for high kmer.  A kmer at this or above may be considered non-error.
eclowthresh=2       (eclt) Threshold for low kmer.  Kmers at this and below are considered errors.
eccmaxqual=127      Do not correct bases with quality above this value.
aec=f               (aggressiveErrorCorrection) Sets more aggressive values of ecr=100, ecclimit=7, echt=16, eclt=3.
cec=f               (conservativeErrorCorrection) Sets more conservative values of ecr=180, ecclimit=2, echt=30, eclt=1, sl=4, pl=4.
meo=f               (markErrorsOnly) Marks errors by reducing quality value of suspected errors; does not correct anything.
mue=t               (markUncorrectableErrors) Marks errors only on uncorrectable reads; requires 'ecc=t'.
overlap=f           (ecco) Error correct by read overlap.

Depth binning parameters:
lowbindepth=10      (lbd) Cutoff for low depth bin.
highbindepth=80     (hbd) Cutoff for high depth bin.
outlow=<file>       Pairs in which both reads have a median below lbd go into this file.
outhigh=<file>      Pairs in which both reads have a median above hbd go into this file.
outmid=<file>       All other pairs go into this file.

Histogram parameters:
hist=<file>         Specify a file to write the input kmer depth histogram.
histout=<file>      Specify a file to write the output kmer depth histogram.
histcol=3           (histogramcolumns) Number of histogram columns, 2 or 3.
pzc=f               (printzerocoverage) Print lines in the histogram with zero coverage.
histlen=1048576     Max kmer depth displayed in histogram.  Also affects statistics displayed, but does not affect normalization.

Peak calling parameters:
peaks=<file>        Write the peaks to this file.  Default is stdout.
minHeight=2         (h) Ignore peaks shorter than this.
minVolume=5         (v) Ignore peaks with less area than this.
minWidth=3          (w) Ignore peaks narrower than this.
minPeak=2           (minp) Ignore peaks with an X-value below this.
maxPeak=BIG         (maxp) Ignore peaks with an X-value above this.
maxPeakCount=8      (maxpc) Print up to this many peaks (prioritizing height).

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx31g"
z2="-Xms31g"
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
	freeRam 31000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

normalize() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.KmerNormalize bits=32 $@"
	echo $CMD >&2
	eval $CMD
}

normalize "$@"
