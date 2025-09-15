#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Jonathan Rood
Last modified October 8, 2024

Description:  Merges paired reads into single reads by overlap detection.
With sufficient coverage, can merge nonoverlapping reads by kmer extension.
Kmer modes (Tadpole or Bloom Filter) require much more memory, and should
be used with the bbmerge-auto.sh script rather than bbmerge.sh.
Please read bbmap/docs/guides/BBMergeGuide.txt for more information.

Usage (interleaved):	bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
Usage (twin files):     bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

Input may be stdin or a file, fasta or fastq, raw or gzipped.

Input parameters:
in=null              Primary input. 'in2' will specify a second file.
interleaved=auto     May be set to true or false to override autodetection of
                     whether the input file as interleaved.
reads=-1             Quit after this many read pairs (-1 means all).

Output parameters:
out=<file>           File for merged reads. 'out2' will specify a second file.
outu=<file>          File for unmerged reads. 'outu2' will specify a second file.
outinsert=<file>     (outi) File to write read names and insert sizes.
outadapter=<file>    (outa) File to write consensus adapter sequences.
outc=<file>          File to write input read kmer cardinality estimate.
ihist=<file>         (hist) Insert length histogram output file.
nzo=t                Only print histogram bins with nonzero values.
showhiststats=t      Print extra header lines with statistical information.
ziplevel=2           Set to 1 (lowest) through 9 (max) to change compression
                     level; lower compression is faster.
ordered=f            Output reads in same order as input.
mix=f                Output both the merged (or mergable) and unmerged reads
                     in the same file (out=).  Useful for ecco mode.

Trimming/Filtering parameters:
qtrim=f              Trim read ends to remove bases with quality below minq.
                     Trims BEFORE merging.
                     Values: t (trim both ends), 
                             f (neither end), 
                             r (right end only), 
                             l (left end only).
qtrim2=f             May be specified instead of qtrim to perform trimming 
                     only if merging is unsuccessful, then retry merging.
trimq=10             Trim quality threshold.  This may be a comma-delimited
                     list (ascending) to try multiple values.
minlength=1          (ml) Reads shorter than this after trimming, but before
                     merging, will be discarded. Pairs will be discarded only
                     if both are shorter.
maxlength=-1         Reads with longer insert sizes will be discarded.
tbo=f                (trimbyoverlap) Trim overlapping reads to remove 
                     rightmost (3') non-overlapping portion, instead of joining.
minavgquality=0      (maq) Reads with average quality below this, after 
                     trimming, will not be attempted to be merged.
maxexpectederrors=0  (mee) If positive, reads with more combined expected 
                     errors than this will not be attempted to be merged.
forcetrimleft=0      (ftl) If nonzero, trim left bases of the read to 
                     this position (exclusive, 0-based).
forcetrimright=0     (ftr) If nonzero, trim right bases of the read 
                     after this position (exclusive, 0-based).
forcetrimright2=0    (ftr2) If positive, trim this many bases on the right end.
forcetrimmod=5       (ftm) If positive, trim length to be equal to 
                     zero modulo this number.
ooi=f                Output only incorrectly merged reads, for testing.
trimpolya=t          Trim trailing poly-A tail from adapter output.  Only 
                     affects outadapter.  This also trims poly-A followed
                     by poly-G, which occurs on NextSeq.

Processing Parameters:
usejni=f             (jni) Do overlapping in C code, which is faster.  Requires
                     compiling the C code; details are in /jni/README.txt.
                     However, the jni path is currently disabled.
merge=t              Create merged reads.  If set to false, you can still 
                     generate an insert histogram.
ecco=f               Error-correct the overlapping part, but don't merge.
trimnonoverlapping=f (tno) Trim all non-overlapping portions, leaving only
                     consensus sequence.  By default, only sequence to the 
                     right of the overlap (adapter sequence) is trimmed.
useoverlap=t         Attempt find the insert size using read overlap.
mininsert=15         Minimum insert size to merge reads.
mininsert0=12        Insert sizes less than this will not be considered.
                     Must be less than or equal to mininsert.
minoverlap=12        Minimum number of overlapping bases to allow merging.
minoverlap0=8        Overlaps shorter than this will not be considered.
                     Must be less than or equal to minoverlap.
minq=9               Ignore bases with quality below this.
maxq=41              Cap output quality scores at this.
quantize=1           Set to a higher number to eliminate some quality scores
                     for a lower output filesize.
entropy=t            Increase the minimum overlap requirement for low-
                     complexity reads.
efilter=6            Ban overlaps with over this many times the expected 
                     number of errors.  Lower is more strict. -1 disables.
pfilter=0.00004      Ban improbable overlaps.  Higher is more strict. 0 will
                     disable the filter; 1 will allow only perfect overlaps.
kfilter=0            Ban overlaps that create kmers with count below
                     this value (0 disables).  If this is used minprob should
                     probably be set to 0.  Requires good coverage.
ouq=f                Calculate best overlap using quality values.
owq=t                Calculate best overlap without using quality values.
usequality=t         If disabled, quality values are completely ignored,
                     both for overlap detection and filtering.  May be useful
                     for data with inaccurate quality values.
iupacton=f           (itn) Change ambiguous IUPAC symbols to N.
adapter=             Specify the adapter sequences used for these reads, if
                     known; this can be a fasta file or a literal sequence.
                     Read 1 and 2 can have adapters specified independently
                     with the adapter1 and adapter2 flags.  adapter=default
                     will use a list of common adapter sequences.

Neural Network Mode Parameters:
nn=t                 Use a neural network for increased merging accuracy.
                     This is highly recommended, but will conflict with
                     strictness and ratiomode flags.  Stringency in nn mode
                     should be adjusted via the cutoff flag instead.
cutoff=0.872857      Merge reads with nn score above this value. Lower will
                     increase the merge rate at the cost of false positives.
net=<file>           Optional network to specify (for developer use); the
                     default is bbmap/resources/bbmerge.bbnet.

Ratio Mode Parameters: 
ratiomode=t          Score overlaps based on the ratio of matching to 
                     mismatching bases.
maxratio=0.09        Max error rate; higher increases merge rate.
ratiomargin=5.5      Lower increases merge rate; min is 1.
ratiooffset=0.55     Lower increases merge rate; min is 0.
maxmismatches=20     Maximum mismatches allowed in overlapping region.
ratiominoverlapreduction=3  This is the difference between minoverlap in 
                     flat mode and minoverlap in ratio mode; generally, 
                     minoverlap should be lower in ratio mode.
minsecondratio=0.1   Cutoff for second-best overlap ratio.
forcemerge=f         Disable all filters and just merge everything 
                     (not recommended).

Strictness (these are mutually exclusive macros that set other parameters):
strict=f             Decrease false positive rate and merging rate.
verystrict=f         (vstrict) Greatly decrease FP and merging rate.
ultrastrict=f        (ustrict) Decrease FP and merging rate even more.
maxstrict=f          (xstrict) Maximally decrease FP and merging rate.
loose=f              Increase false positive rate and merging rate.
veryloose=f          (vloose) Greatly increase FP and merging rate.
ultraloose=f         (uloose) Increase FP and merging rate even more.
maxloose=f           (xloose) Maximally decrease FP and merging rate.
fast=f               Fastest possible mode; less accurate.

Tadpole Parameters (for read extension and error-correction):
*Note: These require more memory and should be run with bbmerge-auto.sh.*
k=31                 Kmer length.  31 (or less) is fastest and uses the least
                     memory, but higher values may be more accurate.  
                     60 tends to work well for 150bp reads.
extend=0             Extend reads to the right this much before merging.
                     Requires sufficient (>5x) kmer coverage.
extend2=0            Extend reads this much only after a failed merge attempt,
                     or in rem/rsem mode.
iterations=1         (ei) Iteratively attempt to extend by extend2 distance
                     and merge up to this many times.
rem=f                (requireextensionmatch) Do not merge if the predicted
                     insert size differs before and after extension.
                     However, if only the extended reads overlap, then that
                     insert will be used.  Requires setting extend2.
rsem=f               (requirestrictextensionmatch) Similar to rem but stricter.
                     Reads will only merge if the predicted insert size before
                     and after extension match.  Requires setting extend2.
                     Enables the lowest possible false-positive rate.
ecctadpole=f         (ecct) If reads fail to merge, error-correct with Tadpole
                     and try again.  This happens prior to extend2.
reassemble=t         If ecct is enabled, use Tadpole's reassemble mode for 
                     error correction.  Alternatives are pincer and tail.
removedeadends       (shave) Remove kmers leading to dead ends.
removebubbles        (rinse) Remove kmers in error bubbles.
mindepthseed=3       (mds) Minimum kmer depth to begin extension.
mindepthextend=2     (mde) Minimum kmer depth continue extension.
branchmult1=20       Min ratio of 1st to 2nd-greatest path depth at high depth.
branchmult2=3        Min ratio of 1st to 2nd-greatest path depth at low depth.
branchlower=3        Max value of 2nd-greatest path depth to be considered low.
ibb=t                Ignore backward branches when extending.
extra=<file>         A file or comma-delimited list of files of reads to use
                     for kmer counting, but not for merging or output.
prealloc=f           Pre-allocate memory rather than dynamically growing; 
                     faster and more memory-efficient for large datasets.  
                     A float fraction (0-1) may be specified, default 1.
prefilter=0          If set to a positive integer, use a countmin sketch to 
                     ignore kmers with depth of that value or lower, to
                     reduce memory usage.
filtermem=0          Allows manually specifying prefilter memory in bytes, for
                     deterministic runs.  0 will set it automatically.
minprob=0.5          Ignore kmers with overall probability of correctness 
                     below this, to reduce memory usage.
minapproxoverlap=26  For rem mode, do not merge reads if the extended reads
                     indicate that the raw reads should have overlapped by
                     at least this much, but no overlap was found.


Bloom Filter Parameters (for kmer operations with less memory than Tadpole)
*Note: These require more memory and should be run with bbmerge-auto.sh.*
eccbloom=f           (eccb) If reads fail to merge, error-correct with bbcms
                     and try again.
testmerge=f          Test kmer counts around the read merge junctions.  If
                     it appears that the merge created new errors, undo it.
                     This reduces the false-positive rate, but not as much as
                     rem or rsem.

Java Parameters:
-Xmx                 This will set Java's memory usage, 
                     overriding autodetection.
                     For example, -Xmx400m will specify 400 MB RAM.
-eoom                This flag will cause the process to exit if an 
                     out-of-memory exception occurs.  Requires Java 8u92+.
-da                  Disable assertions.

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
JNI="-Djava.library.path=""$DIR""jni/"
#JNI=""

z="-Xmx1000m"
z2="-Xms1000m"
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

function merge() {
	local CMD="java $EA $EOOM $z $z2 $JNI $SIMD -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	eval $CMD
}

merge "$@"
