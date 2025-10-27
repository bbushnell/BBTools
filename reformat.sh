#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 22, 2023

Description:  Reformats reads to change ASCII quality encoding, interleaving, file format, or compression format.
Optionally performs additional functions such as quality trimming, subsetting, and subsampling.
Supports fastq, fasta, fasta+qual, scarf, oneline, sam, bam, gzip, bz2.
Please read bbmap/docs/guides/ReformatGuide.txt for more information.

Usage:  reformat.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters and their defaults:

ow=f                    (overwrite) Overwrites files that already exist.
app=f                   (append) Append to files that already exist.
zl=4                    (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f                   (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70            Length of lines in fasta output.
fastareadlen=0          Set to a non-zero number to break fasta files into reads of at most this length.
fastaminlen=1           Ignore fasta reads shorter than this.
qin=auto                ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto               ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
qfake=30                Quality value used for fasta to fastq reformatting.
qfin=<.qual file>       Read qualities from this qual file, for the reads coming from 'in=<fasta file>'
qfin2=<.qual file>      Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'
qfout=<.qual file>      Write qualities from this qual file, for the reads going to 'out=<fasta file>'
qfout2=<.qual file>     Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'
outsingle=<file>        (outs) If a read is longer than minlength and its mate is shorter, the longer one goes here.
deleteinput=f           Delete input upon successful completion.
ref=<file>              Optional reference fasta for sam processing.

Processing Parameters:

verifypaired=f          (vpair) When true, checks reads to see if the names look paired.  Prints an error message if not.
verifyinterleaved=f     (vint) sets 'vpair' to true and 'interleaved' to true.
allowidenticalnames=f   (ain) When verifying pair names, allows identical names, instead of requiring /1 and /2 or 1: and 2:
tossbrokenreads=f       (tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash.
ignorebadquality=f      (ibq) Fix out-of-range quality values instead of crashing with a warning.
addslash=f              Append ' /1' and ' /2' to read names, if not already present.  Please include the flag 'int=t' if the reads are interleaved.
spaceslash=t            Put a space before the slash in addslash mode.
addcolon=f              Append ' 1:' and ' 2:' to read names, if not already present.  Please include the flag 'int=t' if the reads are interleaved.
underscore=f            Change whitespace in read names to underscores.
rcomp=f                 (rc) Reverse-complement reads.
rcompmate=f             (rcm) Reverse-complement read 2 only.
comp=f                  (complement) Reverse-complement reads.
changequality=t         (cq) N bases always get a quality of 0 and ACGT bases get a min quality of 2.
quantize=f              Quantize qualities to a subset of values like NextSeq.  Can also be used with comma-delimited list, like quantize=0,8,13,22,27,32,37
tuc=f                   (touppercase) Change lowercase letters in reads to uppercase.
uniquenames=f           Make duplicate names unique by appending _<number>.
remap=                  A set of pairs: remap=CTGN will transform C>T and G>N.
                        Use remap1 and remap2 to specify read 1 or 2.
iupacToN=f              (itn) Convert non-ACGTN symbols to N.
monitor=f               Kill this process if it crashes.  monitor=600,0.01 would kill after 600 seconds under 1% usage.
crashjunk=t             Crash when encountering reads with invalid bases.
tossjunk=f              Discard reads with invalid characters as bases.
fixjunk=f               Convert invalid bases to N (or X for amino acids).
dotdashxton=f           Specifically convert . - and X to N (or X for amino acids).
recalibrate=f           (recal) Recalibrate quality scores.  Must first generate matrices with CalcTrueQuality.
maxcalledquality=41     Quality scores capped at this upper bound.
mincalledquality=2      Quality scores of ACGT bases will be capped at lower bound.
trimreaddescription=f   (trd) Trim the names of reads after the first whitespace.
trimrname=f             For sam/bam files, trim rname/rnext fields after the first space.
fixheaders=f            Replace characters in headers such as space, *, and | to make them valid file names.
warnifnosequence=t      For fasta, issue a warning if a sequenceless header is encountered.
warnfirsttimeonly=t     Issue a warning for only the first sequenceless header.
utot=f                  Convert U to T (for RNA -> DNA translation).
padleft=0               Pad the left end of sequences with this many symbols.
padright=0              Pad the right end of sequences with this many symbols.
pad=0                   Set padleft and padright to the same value.
padsymbol=N             Symbol to use for padding.

Histogram output parameters:

bhist=<file>            Base composition histogram by position.
qhist=<file>            Quality histogram by position.
qchist=<file>           Count of bases with each quality value.
aqhist=<file>           Histogram of average read quality.
bqhist=<file>           Quality histogram designed for box plots.
lhist=<file>            Read length histogram.
gchist=<file>           Read GC content histogram.
gcbins=100              Number gchist bins.  Set to 'auto' to use read length.
gcplot=f                Add a graphical representation to the gchist.
maxhistlen=6000         Set an upper bound for histogram lengths; higher uses more memory.
                        The default is 6000 for some histograms and 80000 for others.

Histogram parameters for sam files only (requires sam format 1.4 or higher):

ehist=<file>            Errors-per-read histogram.
qahist=<file>           Quality accuracy histogram of error rates versus quality score.
indelhist=<file>        Indel length histogram.
mhist=<file>            Histogram of match, sub, del, and ins rates by read location.
ihist=<file>            Insert size histograms.  Requires paired reads in a sam file.
idhist=<file>           Histogram of read count versus percent identity.
idbins=100              Number idhist bins.  Set to 'auto' to use read length.

Sampling parameters:

reads=-1                Set to a positive number to only process this many INPUT reads (or pairs), then quit.
skipreads=-1            Skip (discard) this many INPUT reads before processing the rest.
samplerate=1            Randomly output only this fraction of reads; 1 means sampling is disabled.
sampleseed=-1           Set to a positive number to use that prng seed for sampling (allowing deterministic sampling).
samplereadstarget=0     (srt) Exact number of OUTPUT reads (or pairs) desired.
samplebasestarget=0     (sbt) Exact number of OUTPUT bases desired.
                        Important: srt/sbt flags should not be used with stdin, samplerate, qtrim, minlength, or minavgquality.
upsample=f              Allow srt/sbt to upsample (duplicate reads) when the target is greater than input.
prioritizelength=f      If true, calculate a length threshold to reach the target, and retain all reads of at least that length (must set srt or sbt).

Trimming and filtering parameters:

qtrim=f                 Trim read ends to remove bases with quality below trimq.
                        Values: t (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window).
trimq=6                 Regions with average quality BELOW this will be trimmed.  Can be a floating-point number like 7.3.
minlength=0             (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
mlf=0                   (mlf) Reads shorter than this fraction of original length after trimming will be discarded.
maxlength=0             If nonzero, reads longer than this after trimming will be discarded.
breaklength=0           If nonzero, reads longer than this will be broken into multiple reads of this length.  Does not work for paired reads.
requirebothbad=t        (rbb) Only discard pairs if both reads are shorter than minlen.
invertfilters=f         (invert) Output failing reads instead of passing reads.
minavgquality=0         (maq) Reads with average quality (after trimming) below this will be discarded.
maqb=0                  If positive, calculate maq from this many initial bases.
chastityfilter=f        (cf) Reads with names  containing ' 1:Y:' or ' 2:Y:' will be discarded.
barcodefilter=f         Remove reads with unexpected barcodes if barcodes is set, or barcodes containing 'N' otherwise.  
                        A barcode must be the last part of the read header.
barcodes=               Comma-delimited list of barcodes or files of barcodes.
maxns=-1                If 0 or greater, reads with more Ns than this (after trimming) will be discarded.
minconsecutivebases=0   (mcb) Discard reads without at least this many consecutive called bases.
forcetrimleft=0         (ftl) If nonzero, trim left bases of the read to this position (exclusive, 0-based).
forcetrimright=-1       (ftr) If nonnegative, trim right bases of the read after this position (exclusive, 0-based).
forcetrimright2=0       (ftr2) If positive, trim this many bases on the right end.
forcetrimmod=5          (ftm) If positive, trim length to be equal to zero modulo this number.
mingc=0                 Discard reads with GC content below this.
maxgc=1                 Discard reads with GC content above this.
gcpairs=t               Use average GC of paired reads.
                        Also affects gchist.

Tag-filtering parameters:

tag=                    Look for this tag in the header to filter by the next value.  To filter reads
                        with a header like 'foo,depth=5.5,bar' where you only want depths
                        of at least 3, the necessary flags would be 'tag=depth= minvalue=3 delimiter=,'
delimiter=              Character after the end of the value, such as delimiter=X.  Control and
                        whitespace symbols may be spelled out, like delimiter=tab or delimiter=pipe.
                        The tag may contain the delimiter.  If the value is the last term in the header,
                        the delimiter doesn't matter but is still required.
minvalue=               If set, only accept a numeric value of at least this.
maxvalue=               If set, only accept a numeric value of at most this.
value=                  If set, only accept a string value of exactly this.

Illumina-specific parameters:
top=true                Include reads from the top of the flowcell.
bottom=true             Include reads from the bottom of the flowcell.

Sam and bam processing parameters:

mappedonly=f            Toss unmapped reads.
unmappedonly=f          Toss mapped reads.
pairedonly=f            Toss reads that are not mapped as proper pairs.
unpairedonly=f          Toss reads that are mapped as proper pairs.
primaryonly=f           Toss secondary alignments.  Set this to true for sam to fastq conversion.
minmapq=-1              If non-negative, toss reads with mapq under this.
maxmapq=-1              If non-negative, toss reads with mapq over this.
requiredbits=0          (rbits) Toss sam lines with any of these flag bits unset.  Similar to samtools -f.
filterbits=0            (fbits) Toss sam lines with any of these flag bits set.  Similar to samtools -F.
stoptag=f               Set to true to write a tag indicating read stop location, prefixed by YS:i:
sam=                    Set to 'sam=1.3' to convert '=' and 'X' cigar symbols (from sam 1.4+ format) to 'M'.
                        Set to 'sam=1.4' to convert 'M' to '=' and 'X' (sam=1.4 requires MD tags to be present, or ref to be specified).

Sam and bam alignment filtering parameters:
These require = and X symbols in cigar strings, or MD tags, or a reference fasta.
-1 means disabled; to filter reads with any of a symbol type, set to 0.

subfilter=-1            Discard reads with more than this many substitutions.
minsubs=-1              Discard reads with fewer than this many substitutions.
insfilter=-1            Discard reads with more than this many insertions.
delfilter=-1            Discard reads with more than this many deletions.
indelfilter=-1          Discard reads with more than this many indels.
editfilter=-1           Discard reads with more than this many edits.
inslenfilter=-1         Discard reads with an insertion longer than this.
dellenfilter=-1         Discard reads with a deletion longer than this.
minidfilter=-1.0        Discard reads with identity below this (0-1).
maxidfilter=1.0         Discard reads with identity above this (0-1).
clipfilter=-1           Discard reads with more than this many soft-clipped bases.

Kmer counting and cardinality estimation parameters:
k=0                     If positive, count the total number of kmers.
cardinality=f           (loglog) Count unique kmers using the LogLog algorithm.
loglogbuckets=1999      Use this many buckets for cardinality estimation.

Shortcuts: 
The # symbol will be substituted for 1 and 2.  The % symbol in out will be substituted for input name minus extensions.
For example:
reformat.sh in=read#.fq out=%.fa
...is equivalent to:
reformat.sh in1=read1.fq in2=read2.fq out1=read1.fa out2=read2.fa

Java Parameters:
-Xmx                    This will set Java's memory usage, overriding autodetection.
                        -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                        The max is typically 85% of physical memory.
-eoom                   This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da                     Disable assertions.

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

z="-Xmx300m"
z2="-Xms300m"
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

function reformat() {
	local CMD="java $EA $EOOM $z $z2 $SIMD -cp $CP jgi.ReformatReads $@"
	echo $CMD >&2
	eval $CMD
}

reformat "$@"
