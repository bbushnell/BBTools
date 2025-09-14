#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 10, 2025

Description:  Sorts sequences to put similar reads near each other.
Can be used for increased compression or error correction.
Please read bbmap/docs/guides/ClumpifyGuide.txt for more information.

Usage:   clumpify.sh in=<file> out=<file> reorder

Input may be fasta or fastq, compressed or uncompressed.  Cannot accept sam.

Parameters and their defaults:
in=<file>           Input file.
in2=<file>          Optional input for read 2 of twin paired files.
out=<file>          Output file.  May not be standard out.
out2=<file>         Optional output for read 2 of twin paired files.
groups=auto         Use this many intermediate files (to save memory).
                    1 group is fastest.  Auto will estimate the number
                    of groups needed based on the file size, so it should
                    not ever run out of memory.
lowcomplexity=f     For compressed low-complexity libraries such as RNA-seq,
                    this will more conservatively estimate how much memory
                    is needed to automatically decide the number of groups.
rcomp=f             Give read clumps the same orientation to increase 
                    compression.  Should be disabled for paired reads.
overwrite=f         (ow) Set to false to force the program to abort rather 
                    than overwrite an existing file.
qin=auto            Auto-detect input quality encoding.  May be set to:
                       33:  ASCII-33 (Sanger) encoding.
                       64:  ASCII-64 (old Illumina) encoding.
                    All modern sequence is encoded as ASCII-33.
qout=auto           Use input quality encoding as output quality encoding.
changequality=f     (cq) If true, fix broken quality scores such as Ns with
                    Q>0.  Default is false to ensure lossless compression.
fastawrap=70        Set to a higher number like 4000 for longer lines in 
                    fasta format, which increases compression.

Compression parameters:
ziplevel=6          (zl) Gzip compression level (1-11).  Higher is slower.
                    Level 11 is only available if pigz is installed and is
                    extremely slow to compress, but faster to decompress.
                    Naming the output file to *.bz2 will use bzip2 instead of
                    gzip for ~9% additional compression, which requires
                    bzip2, pbzip2, or lbzip2 in the path.
blocksize=128       Size of blocks for pigz, in kb.  Higher gives slightly
                    better compression.
shortname=f         Make the names as short as possible.  'shortname=shrink'
                    will shorten the names where possible, but retain the 
                    flowcell and barcode information.
reorder=f           Reorder clumps for additional compression.  Only valid
                    when groups=1, passes=1, and ecc=f.  Possible modes:
                       f:  Do not reorder clumps.
                       c:  Reorder using consensus reads.  Uses additional
                           time and memory.
                       p:  Reorder using pair information.  Requires paired
                           reads.  Yields the highest compression.
                       a:  Automatically choose between 'c' and 'p'.  The
                           flag reorder with no argument will set 'reorder=a'.
quantize=f          Bin the quality scores, like NextSeq.  This greatly
                    increases compression, but information is lost.

Temp file parameters:
compresstemp=auto   (ct) Gzip temporary files.  By default temp files will be
                    compressed if the output file is compressed.
deletetemp=t        Delete temporary files.
deleteinput=f       Delete input upon successful completion.
usetmpdir=f         Use tmpdir for temp files.
tmpdir=             By default, this is the environment variable TMPDIR.

Hashing parameters:
k=31                Use kmers of this length (1-31).  Shorter kmers may
                    increase compression, but 31 is recommended for error
                    correction.
mincount=0          Don't use pivot kmers with count less than this.
                    Setting mincount=2 can increase compression.
                    Increases time and memory usage.
seed=1              Random number generator seed for hashing.  
                    Set to a negative number to use a random seed.
hashes=4            Use this many masks when hashing.  0 uses raw kmers.
                    Often hashes=0 increases compression, but it should
                    not be used with error-correction.
border=1            Do not use kmers within this many bases of read ends.

Deduplication parameters:
dedupe=f            Remove duplicate reads.  For pairs, both must match.
                    By default, deduplication does not occur.
                    If dedupe and markduplicates are both false, none of
                    the other duplicate-related flags will have any effect.
markduplicates=f    Don't remove; just append ' duplicate' to the name.
allduplicates=f     Mark or remove all copies of duplicates, instead of
                    keeping the highest-quality copy.
addcount=f          Append the number of copies to the read name.
                    Mutually exclusive with markduplicates or allduplicates.
entryfilter=f       This assists in removing exact duplicates, which saves
                    memory in libraries that split unevenly due to huge
                    numbers of duplicates.  Enabled automatically as needed.
subs=2              (s) Maximum substitutions allowed between duplicates.
subrate=0.0         (dsr) If set, the number of substitutions allowed will be
                    max(subs, subrate*min(length1, length2)) for 2 sequences.
allowns=t           No-called bases will not be considered substitutions.
scanlimit=5         (scan) Continue for this many reads after encountering a
                    nonduplicate.  Improves detection of inexact duplicates.
umi=f               If reads have UMIs in the headers, require them to match
                    to consider the reads duplicates.
umisubs=0           Consider UMIs as matching if they have up to this many
                    mismatches.
containment=f       Allow containments (where one sequence is shorter).
affix=f             For containments, require one sequence to be an affix
                    (prefix or suffix) of the other.
optical=f           If true, mark or remove optical duplicates only.
                    This means they are Illumina reads within a certain
                    distance on the flowcell.  Normal Illumina names needed.
                    Also for tile-edge and well duplicates.
dupedist=40         (dist) Max distance to consider for optical duplicates.
                    Higher removes more duplicates but is more likely to
                    remove PCR rather than optical duplicates.
                    This is platform-specific; recommendations:
                       NextSeq      40  (and spany=t)
                       HiSeq 1T     40
                       HiSeq 2500   40
                       HiSeq 3k/4k  2500
                       Novaseq6000  12000
                       NovaseqX+    50
spany=f             Allow reads to be considered optical duplicates if they
                    are on different tiles, but are within dupedist in the
                    y-axis.  Should only be enabled when looking for 
                    tile-edge duplicates (as in NextSeq).
spanx=f             Like spany, but for the x-axis.  Not necessary 
                    for NextSeq.
spantiles=f         Set both spanx and spany.
adjacent=f          Limit tile-spanning to adjacent tiles (those with 
                    consecutive numbers).
*** Thus, for NextSeq, the recommended deduplication flags are: ***
dedupe optical spany adjacent

Pairing/ordering parameters (for use with error-correction):
unpair=f            For paired reads, clump all of them rather than just
                    read 1.  Destroys pairing.  Without this flag, for paired
                    reads, only read 1 will be error-corrected.
repair=f            After clumping and error-correction, restore pairing.
                    If groups>1 this will sort by name which will destroy
                    clump ordering; with a single group, clumping will
                    be retained.

Error-correction parameters:
ecc=f               Error-correct reads.  Requires multiple passes for
                    complete correction.
ecco=f              Error-correct paired reads via overlap before clumping.
passes=1            Use this many error-correction passes.  6 passes are 
                    suggested, though more will be more through.
conservative=f      Only correct the highest-confidence errors, to minimize
                    chances of eliminating a minor allele or inexact repeat.
aggressive=f        Maximize the number of errors corrected.
consensus=f         Output consensus sequence instead of clumps.

Advanced error-correction parameters:
mincc=4             (mincountcorrect) Do not correct to alleles occuring less
                    often than this.
minss=4             (minsizesplit) Do not split into new clumps smaller than 
                    this.
minsfs=0.17         (minsizefractionsplit) Do not split on pivot alleles in
                    areas with local depth less than this fraction of clump size.
minsfc=0.20         (minsizefractioncorrect) Do not correct in areas with local
                    depth less than this.
minr=30.0           (minratio) Correct to the consensus if the ratio of the
                    consensus allele to second-most-common allele is >=minr,
                    for high depth.  Actual ratio used is:
                    min(minr, minro+minorCount*minrm+quality*minrqm).
minro=1.9           (minratiooffset) Base ratio.
minrm=1.8           (minratiomult) Ratio multiplier for secondary allele count.
minrqm=0.08         (minratioqmult) Ratio multiplier for base quality.
minqr=2.8           (minqratio) Do not correct bases when cq*minqr>rqsum.
minaqr=0.70         (minaqratio) Do not correct bases when cq*minaqr>5+rqavg.
minid=0.97          (minidentity) Do not correct reads with identity to 
                    consensus less than this.
maxqadjust=0        Adjust quality scores by at most maxqadjust per pass.
maxqi=-1            (maxqualityincorrect) Do not correct bases with quality 
                    above this (if positive).
maxci=-1            (maxcountincorrect) Do not correct alleles with count 
                    above this (if positive).
findcorrelations=t  Look for correlated SNPs in clumps to split into alleles.
maxcorrelations=12  Maximum number of eligible SNPs per clump to consider for
                    correlations.  Increasing this number can reduce false-
                    positive corrections at the possible expense of speed.

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

z="-Xmx2g"
z2="-Xms2g"
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
	#Note that this uses slighly less (82%) of memory than normal to account for multiple pigz instances.
	freeRam 2000m 82
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

clumpify() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP clump.Clumpify $@"
	echo $CMD >&2
	eval $CMD
}

clumpify "$@"
