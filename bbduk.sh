#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 29, 2025

Description:  Compares reads to the kmers in a reference dataset, optionally 
allowing an edit distance. Splits the reads into two outputs - those that 
match the reference, and those that don't. Can also trim (remove) the matching 
parts of the reads rather than binning the reads.
Please read bbmap/docs/guides/BBDukGuide.txt for more information.

Usage:  bbduk.sh in=<input file> out=<output file> ref=<contaminant files>

Input may be stdin or a fasta or fastq file, compressed or uncompressed.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped 
fasta input, set in=stdin.fa.gz

Input parameters:
in=<file>           Main input. in=stdin.fq will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
ref=<file,file>     Comma-delimited list of reference files.
                    In addition to filenames, you may also use the keywords:
                    adapters, artifacts, phix, lambda, pjet, mtst, kapa
literal=<seq,seq>   Comma-delimited list of literal reference sequences.
                    Polymers are also allowed with the 'poly' prefix;
                    for example, 'literal=ATGGT,polyGC' will add both ATGGT
                    and GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC - 32+ of them,
                    enough replicates to ensure that all kmers are present.
touppercase=f       (tuc) Change all bases upper-case.
interleaved=auto    (int) t/f overrides interleaved autodetection.
                    Must be set mainually when streaming fastq input.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
reads=-1            If positive, quit after processing X reads or pairs.
copyundefined=f     (cu) Process non-AGCT IUPAC reference bases by making all
                    possible unambiguous copies.  Intended for short motifs
                    or adapter barcodes, as time/memory use is exponential.
samplerate=1        Set lower to only process a fraction of input reads.
samref=<file>       Optional reference fasta for processing sam files.

Output parameters:
out=<file>          (outnonmatch) Write reads here that do not contain 
                    kmers matching the database.  'out=stdout.fq' will pipe 
                    to standard out.
out2=<file>         (outnonmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outm=<file>         (outmatch) Write reads here that fail filters.  In default
                    kfilter mode, this means any read with a matching kmer.
                    In any mode, it also includes reads that fail filters such
                    as minlength, mingc, maxgc, entropy, etc.  In other words,
                    it includes all reads that do not go to 'out'.
outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outs=<file>         (outsingle) Use this to write singleton reads whose mate 
                    was trimmed shorter than minlen.
stats=<file>        Write statistics about which contamininants were detected.
refstats=<file>     Write statistics on a per-reference-file basis.
rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
dump=<file>         Dump kmer tables to a file, in fasta format.
duk=<file>          Write statistics in duk's format. *DEPRECATED*
nzo=t               Only write statistics about ref sequences with nonzero hits.
overwrite=t         (ow) Grant permission to overwrite files.
showspeed=t         (ss) 'f' suppresses display of processing speed.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=70        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
statscolumns=3      (cols) Number of columns for stats output, 3 or 5.
                    5 includes base counts.
rename=f            Rename reads to indicate which sequences they matched.
refnames=f          Use names of reference files rather than scaffold IDs.
trd=f               Truncate read and ref names at the first whitespace.
ordered=f           Set to true to output reads in same order as input.
maxbasesout=-1      If positive, quit after writing approximately this many
                    bases to out (outu/outnonmatch).
maxbasesoutm=-1     If positive, quit after writing approximately this many
                    bases to outm (outmatch).
json=f              Print to screen in json format.

Histogram output parameters:
bhist=<file>        Base composition histogram by position.
qhist=<file>        Quality histogram by position.
qchist=<file>       Count of bases with each quality value.
aqhist=<file>       Histogram of average read quality.
bqhist=<file>       Quality histogram designed for box plots.
lhist=<file>        Read length histogram.
phist=<file>        Polymer length histogram.
gchist=<file>       Read GC content histogram.
enthist=<file>      Read entropy histogram.
ihist=<file>        Insert size histogram, for paired reads in mapped sam.
gcbins=100          Number gchist bins.  Set to 'auto' to use read length.
maxhistlen=6000     Set an upper bound for histogram lengths; higher uses 
                    more memory.  The default is 6000 for some histograms
                    and 80000 for others.

Histogram parameters for mapped sam/bam files only:
histbefore=t        Calculate histograms from reads before processing.
ehist=<file>        Errors-per-read histogram.
qahist=<file>       Quality accuracy histogram of error rates versus quality 
                    score.
indelhist=<file>    Indel length histogram.
mhist=<file>        Histogram of match, sub, del, and ins rates by position.
idhist=<file>       Histogram of read count versus percent identity.
idbins=100          Number idhist bins.  Set to 'auto' to use read length.
varfile=<file>      Ignore substitution errors listed in this file when 
                    calculating error rates.  Can be generated with
                    CallVariants.
vcf=<file>          Ignore substitution errors listed in this VCF file 
                    when calculating error rates.
ignorevcfindels=t   Also ignore indels listed in the VCF.

Processing parameters:
k=31                Kmer length used for finding contaminants.  Contaminants 
                    shorter than k will not be found.  k must be at least 1.
rcomp=t             Look for reverse-complements of kmers in addition to 
                    forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to 
                    increase sensitivity in the presence of errors.  This may
                    also be set to a number, e.g. mm=3, to mask that many bp.
                    The default mm=t corresponds to mm=1 for odd-length kmers
                    and mm=2 for even-length kmers (as of v39.04), while
                    mm=f is always equivalent to mm=0.
minkmerhits=1       (mkh) Reads need at least this many matching kmers 
                    to be considered as matching the reference.
minkmerfraction=0.0 (mkf) A reads needs at least this fraction of its total
                    kmers to hit a ref, in order to be considered a match.
                    If this and minkmerhits are set, the greater is used.
mincovfraction=0.0  (mcf) A reads needs at least this fraction of its total
                    bases to be covered by ref kmers to be considered a match.
                    If specified, mcf overrides mkh and mkf.
hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
                    Memory use is proportional to (3*K)^hdist.
qhdist=0            Hamming distance for query kmers; impacts speed, not memory.
editdistance=0      (edist) Maximum edit distance from ref kmers (subs 
                    and indels).  Memory use is proportional to (8*K)^edist.
hammingdistance2=0  (hdist2) Sets hdist for short kmers, when using mink.
qhdist2=0           Sets qhdist for short kmers, when using mink.
editdistance2=0     (edist2) Sets edist for short kmers, when using mink.
forbidn=f           (fn) Forbids matching of read kmers containing N.
                    By default, these will match a reference 'A' if 
                    hdist>0 or edist>0, to increase sensitivity.
removeifeitherbad=t (rieb) Paired reads get sent to 'outmatch' if either is 
                    match (or either is trimmed shorter than minlen).  
                    Set to false to require both.
trimfailures=f      Instead of discarding failed reads, trim them to 1bp.
                    This makes the statistics a bit odd.
findbestmatch=f     (fbm) If multiple matches, associate read with sequence 
                    sharing most kmers.  Reduces speed.
skipr1=f            Don't do kmer-based operations on read 1.
skipr2=f            Don't do kmer-based operations on read 2.
ecco=f              For overlapping paired reads only.  Performs error-
                    correction with BBMerge prior to kmer operations.
recalibrate=f       (recal) Recalibrate quality scores.  Requires calibration
                    matrices generated by CalcTrueQuality.
sam=<file,file>     If recalibration is desired, and matrices have not already
                    been generated, BBDuk will create them from the sam file.
amino=f             Run in amino acid mode.  Some features have not been
                    tested, but kmer-matching works fine.  Maximum k is 12.

Speed and Memory parameters:
threads=auto        (t) Set number of threads to use; default is number of 
                    logical processors.
prealloc=f          Preallocate memory in table.  Allows faster table loading 
                    and more efficient memory usage, for a large reference.
monitor=f           Kill this process if it crashes.  monitor=600,0.01 would 
                    kill after 600 seconds under 1% usage.
minrskip=1          (mns) Force minimal skip interval when indexing reference 
                    kmers.  1 means use all, 2 means use every other kmer, etc.
maxrskip=1          (mxs) Restrict maximal skip interval when indexing 
                    reference kmers. Normally all are used for scaffolds<100kb, 
                    but with longer scaffolds, up to maxrskip-1 are skipped.
rskip=              Set both minrskip and maxrskip to the same value.
                    If not set, rskip will vary based on sequence length.
qskip=1             Skip query kmers to increase speed.  1 means use all.
speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                    reads and reference.  Increases speed and reduces memory.
Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

Trimming/Filtering/Masking parameters:
Note - if ktrim, kmask, and ksplit are unset, the default behavior is kfilter.
All kmer processing modes are mutually exclusive.
Reads only get sent to 'outm' purely based on kmer matches in kfilter mode.

ktrim=f             Trim reads to remove bases matching reference kmers, plus
                    all bases to the left or right.
                    Values:
                       f (don't trim), 
                       r (trim to the right), 
                       l (trim to the left)
ktrimtips=0         Set this to a positive number to perform ktrim on both
                    ends, examining only the outermost X bases.
kmask=              Replace bases matching ref kmers with another symbol.
                    Allows any non-whitespace character, and processes short
                    kmers on both ends if mink is set.  'kmask=lc' will
                    convert masked bases to lowercase.
maskfullycovered=f  (mfc) Only mask bases that are fully covered by kmers.
ksplit=f            For single-ended reads only.  Reads will be split into
                    pairs around the kmer.  If the kmer is at the end of the
                    read, it will be trimmed instead.  Singletons will go to
                    out, and pairs will go to outm.  Do not use ksplit with
                    other operations such as quality-trimming or filtering.
mink=0              Look for shorter kmers at read tips down to this length, 
                    when k-trimming or masking.  0 means disabled.  Enabling
                    this will disable maskmiddle.
qtrim=f             Trim read ends to remove bases with quality below trimq.
                    Performed AFTER looking for kmers.  Values: 
                       rl (trim both ends), 
                       f (neither end), 
                       r (right end only), 
                       l (left end only),
                       w (sliding window).
trimq=6             Regions with average quality BELOW this will be trimmed,
                    if qtrim is set to something other than f.  Can be a 
                    floating-point number like 7.3.
quantize            Bin quality scores to reduce file size.  quantize=2 will
                    eliminate all odd quality scores, while quantize=0,10,37
                    will only allow qualty scores of 0, 10, or 37.
trimclip=f          Trim soft-clipped bases from sam files.
minlength=10        (ml) Reads shorter than this after trimming will be 
                    discarded.  Pairs will be discarded if both are shorter.
mlf=0               (minlengthfraction) Reads shorter than this fraction of 
                    original length after trimming will be discarded.
maxlength=          Reads longer than this after trimming will be discarded.
minavgquality=0     (maq) Reads with average quality (after trimming) below 
                    this will be discarded.
maqb=0              If positive, calculate maq from this many initial bases.
minbasequality=0    (mbq) Reads with any base below this quality (after 
                    trimming) will be discarded.
maxns=-1            If non-negative, reads with more Ns than this 
                    (after trimming) will be discarded.
mcb=0               (minconsecutivebases) Discard reads without at least 
                    this many consecutive called bases.
ottm=f              (outputtrimmedtomatch) Output reads trimmed to shorter 
                    than minlength to outm rather than discarding.
tp=0                (trimpad) Trim this much extra around matching kmers.
tbo=f               (trimbyoverlap) Trim adapters based on where paired 
                    reads overlap.
strictoverlap=t     Adjust sensitivity for trimbyoverlap mode.
minoverlap=14       Require this many bases of overlap for detection.
mininsert=40        Require insert size of at least this for overlap.
                    Should be reduced to 16 for small RNA sequencing.
tpe=f               (trimpairsevenly) When kmer right-trimming, trim both 
                    reads to the minimum length of either.
forcetrimleft=0     (ftl) If positive, trim bases to the left of this position
                    (exclusive, 0-based).
forcetrimright=0    (ftr) If positive, trim bases to the right of this position
                    (exclusive, 0-based).
forcetrimright2=0   (ftr2) If positive, trim this many bases on the right end.
forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
                    modulo this number.
restrictleft=0      If positive, only look for kmer matches in the 
                    leftmost X bases.
restrictright=0     If positive, only look for kmer matches in the 
                    rightmost X bases.
NOTE:  restrictleft and restrictright are mutually exclusive.  If trimming
       both ends is desired, use ktrimtips.
mingc=0             Discard reads with GC content below this.
maxgc=1             Discard reads with GC content above this.
gcpairs=t           Use average GC of paired reads.
                    Also affects gchist.
tossjunk=f          Discard reads with invalid characters as bases.
swift=f             Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.

Header-parsing parameters - these require Illumina headers:
chastityfilter=f    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
barcodefilter=f     Remove reads with unexpected barcodes if barcodes is set,
                    or barcodes containing 'N' otherwise.  A barcode must be
                    the last part of the read header.  Values:
                       t:     Remove reads with bad barcodes.
                       f:     Ignore barcodes.
                       crash: Crash upon encountering bad barcodes.
barcodes=           Comma-delimited list of barcodes or files of barcodes.
xmin=-1             If positive, discard reads with a lesser X coordinate.
ymin=-1             If positive, discard reads with a lesser Y coordinate.
xmax=-1             If positive, discard reads with a greater X coordinate.
ymax=-1             If positive, discard reads with a greater Y coordinate.

Polymer trimming parameters:
trimpolya=0         If greater than 0, trim poly-A or poly-T tails of
                    at least this length on either end of reads.
trimpolygleft=0     If greater than 0, trim poly-G prefixes of at least this
                    length on the left end of reads.  Does not trim poly-C.
trimpolygright=0    If greater than 0, trim poly-G tails of at least this 
                    length on the right end of reads.  Does not trim poly-C.
trimpolyg=0         This sets both left and right at once.
filterpolyg=0       If greater than 0, remove reads with a poly-G prefix of
                    at least this length (on the left).
Note: there are also equivalent poly-C flags.

Polymer tracking parameters:
pratio=base,base    'pratio=G,C' will print the ratio of G to C polymers.
plen=20             Length of homopolymers to count.

Entropy/Complexity parameters:
entropy=-1          Set between 0 and 1 to filter reads with entropy below
                    that value.  Higher is more stringent.
entropywindow=50    Calculate entropy using a sliding window of this length.
entropyk=5          Calculate entropy using kmers of this length.
minbasefrequency=0  Discard reads with a minimum base frequency below this.
entropytrim=f       Values:
                       f:  (false) Do not entropy-trim.
                       r:  (right) Trim low entropy on the right end only.
                       l:  (left) Trim low entropy on the left end only.
                       rl: (both) Trim low entropy on both ends.
entropymask=f       Values:
                       f:  (filter) Discard low-entropy sequences.
                       t:  (true) Mask low-entropy parts of sequences with N.
                       lc: Change low-entropy parts of sequences to lowercase.
entropymark=f       Mark each base with its entropy value.  This is on a scale
                    of 0-41 and is reported as quality scores, so the output
                    should be fastq or fasta+qual.
NOTE: If set, entropytrim overrides entropymask.

Cardinality estimation parameters:
cardinality=f       (loglog) Count unique kmers using the LogLog algorithm.
cardinalityout=f    (loglogout) Count unique kmers in output reads.
loglogk=31          Use this kmer length for counting.
loglogbuckets=2048  Use this many buckets for counting.
khist=<file>        Kmer frequency histogram; plots number of kmers versus
                    kmer depth.  This is approximate.
khistout=<file>     Kmer frequency histogram for output reads.

Side Channel Parameters:
sideout=<file>      Output for aligned reads.
sideref=phix        Reference for side-channel alignment; must be a single
                    sequence and virtually repeat-free at selected k.
sidek1=17           Kmer length for seeding alignment to reference.
sidek2=13           Kmer length for seeding alignment of unaligned reads
                    with an aligned mate.
sideminid1=0.66     Minimum identity to accept individual alignments.
sideminid2=0.58     Minimum identity for aligning reads with aligned mates.
sidemm1=1           Middle mask length for sidek1.
sidemm2=1           Middle mask length for sidek2.
Note:  The side channel is a special additional output that allows alignment
to a secondary reference while also doing trimming.  Alignment does not affect
whether reads go to the normal outputs (out, outm).  The main purpose is to
simplify pipelines that need trimmed, aligned phiX reads for recalibration.


Java Parameters:

-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

z="-Xmx1g"
z2="-Xms1g"
set=0
silent=0
json=0

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
	freeRam 1400m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbduk() {
	local CMD="java $EA $EOOM $z $z2 $SIMD $JNI -cp $CP jgi.BBDuk $@"
	if [[ $silent == 0 ]] && [[ $json == 0 ]]; then
		echo $CMD >&2
	fi
	eval $CMD
}

bbduk "$@"
