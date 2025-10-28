#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 3, 2021

Description:  Uses kmer counts to assemble contigs, extend sequences, 
or error-correct reads.  Tadpole has no upper bound for kmer length,
but some values are not supported.  Specifically, it allows 1-31,
multiples of 2 from 32-62, multiples of 3 from 63-93, etc.
Please read bbmap/docs/guides/TadpoleGuide.txt for more information.

Usage (Assembly):  tadpole.sh k=62 in=<reads> out=<contigs>
Extension:    tadpole.sh k=62 in=<reads> out=<extended> mode=extend
Correction:   tadpole.sh k=62 in=<reads> out=<corrected> mode=correct

Recommended parameters for optimal assembly:
tadpole.sh in=<reads> out=<contigs> shave rinse pop k=<50-70% of read length>

Extension and correction may be done simultaneously.  Error correction on 
multiple files may be done like this:

tadpole.sh in=libA_r1.fq,libA_merged.fq in2=libA_r2.fq,null extra=libB_r1.fq out=ecc_libA_r1.fq,ecc_libA_merged.fq out2=ecc_libA_r2.fq,null mode=correct

Extending contigs with reads could be done like this:

tadpole.sh in=contigs.fa out=extended.fa el=100 er=100 mode=extend extra=reads.fq k=62


Input parameters:
in=<file>           Primary input file for reads to use as kmer data.
in2=<file>          Second input file for paired data.
extra=<file>        Extra files for use as kmer data, but not for error-
                    correction or extension.
reads=-1            Only process this number of reads, then quit (-1 means all).
NOTE: in, in2, and extra may also be comma-delimited lists of files.

Output parameters:
out=<file>          Write contigs (in contig mode) or corrected/extended 
                    reads (in other modes).
out2=<file>         Second output file for paired output.
outd=<file>         Write discarded reads, if using junk-removal flags.
dot=<file>          Write a contigs connectivity graph (partially implemented)
dump=<file>         Write kmers and their counts.
fastadump=t         Write kmers and counts as fasta versus 2-column tsv.
mincounttodump=1    Only dump kmers with at least this depth.
showstats=t         Print assembly statistics after writing contigs.

Prefiltering parameters:
prefilter=0         If set to a positive integer, use a countmin sketch
                    to ignore kmers with depth of that value or lower.
prehashes=2         Number of hashes for prefilter.
prefiltersize=0.2   (pff) Fraction of memory to use for prefilter.
minprobprefilter=t  (mpp) Use minprob for the prefilter.
prepasses=1         Use this many prefiltering passes; higher be more thorough
                    if the filter is very full.  Set to 'auto' to iteratively 
                    prefilter until the remaining kmers will fit in memory.
onepass=f           If true, prefilter will be generated in same pass as kmer
                    counts.  Much faster but counts will be lower, by up to
                    prefilter's depth limit.
filtermem=0         Allows manually specifying prefilter memory in bytes, for
                    deterministic runs.  0 will set it automatically.

Hashing parameters:
k=31                Kmer length (1 to infinity).  Memory use increases with K.
prealloc=t          Pre-allocate memory rather than dynamically growing; 
                    faster and more memory-efficient.  A float fraction (0-1)
                    may be specified; default is 1.
minprob=0.5         Ignore kmers with overall probability of correctness below this.
minprobmain=t       (mpm) Use minprob for the primary kmer counts.
threads=X           Spawn X worker threads; default is number of logical processors.
buildthreads=X      Spawn X contig-building threads. If not set, defaults to the same
                    as threads.  Setting this to 1 will make contigs deterministic.
rcomp=t             Store and count each kmer together and its reverse-complement.
coremask=t          All kmer extensions share the same hashcode.
fillfast=t          Speed up kmer extension lookups.

Assembly parameters:
mincountseed=3      (mcs) Minimum kmer count to seed a new contig or begin extension.
mincountextend=2    (mce) Minimum kmer count continue extension of a read or contig.
                    It is recommended that mce=1 for low-depth metagenomes.
mincountretain=0    (mincr) Discard kmers with count below this.
maxcountretain=INF  (maxcr) Discard kmers with count above this.
branchmult1=20      (bm1) Min ratio of 1st to 2nd-greatest path depth at high depth.
branchmult2=3       (bm2) Min ratio of 1st to 2nd-greatest path depth at low depth.
branchlower=3       (blc) Max value of 2nd-greatest path depth to be considered low.
minextension=2      (mine) Do not keep contigs that did not extend at least this much.
mincontig=auto      (minc) Do not write contigs shorter than this.
mincoverage=1       (mincov) Do not write contigs with average coverage below this.
maxcoverage=inf     (maxcov) Do not write contigs with average coverage above this.
trimends=0          (trim) Trim contig ends by this much.  Trimming by K/2 
                    may yield more accurate genome size estimation.
trimcircular=t      Trim one end of contigs ending in LOOP/LOOP by K-1,
                    to eliminate the overlapping portion.
contigpasses=16     Build contigs with decreasing seed depth for this many iterations.
contigpassmult=1.7  Ratio between seed depth of two iterations.
ownership=auto      For concurrency; do not touch.
processcontigs=f    Explore the contig connectivity graph.
popbubbles=t        (pop) Pop bubbles; increases contiguity.  Requires 
                    additional time and memory and forces processcontigs=t.

Processing mode parameters:
mode=contig         contig: Make contigs from kmers.
                    extend: Extend sequences to be longer, and optionally
                            perform error correction.
                    correct: Error correct only.
                    insert: Measure insert sizes.
                    discard: Discard low-depth reads, without error correction.

Extension parameters:
extendleft=100      (el) Extend to the left by at most this many bases.
extendright=100     (er) Extend to the right by at most this many bases.
ibb=t               (ignorebackbranches) Do not stop at backward branches.
extendrollback=3    Trim a random number of bases, up to this many, on reads
                    that extend only partially.  This prevents the creation
                    of sharp coverage discontinuities at branches.

Error-correction parameters:
ecc=f               Error correct via kmer counts.
reassemble=t        If ecc is enabled, use the reassemble algorithm.
pincer=f            If ecc is enabled, use the pincer algorithm.
tail=f              If ecc is enabled, use the tail algorithm.
eccfull=f           If ecc is enabled, use tail over the entire read.
aggressive=f        (aecc) Use aggressive error correction settings.
                    Overrides some other flags like errormult1 and deadzone.
conservative=f      (cecc) Use conservative error correction settings.
                    Overrides some other flags like errormult1 and deadzone.
rollback=t          Undo changes to reads that have lower coverage for
                    any kmer after correction.
markbadbases=0      (mbb) Any base fully covered by kmers with count below 
                    this will have its quality reduced.
markdeltaonly=t     (mdo) Only mark bad bases adjacent to good bases.
meo=t               (markerrorreadsonly) Only mark bad bases in reads 
                    containing errors.
markquality=0       (mq) Set quality scores for marked bases to this.
                    A level of 0 will also convert the base to an N.
errormult1=16       (em1) Min ratio between kmer depths to call an error.
errormult2=2.6      (em2) Alternate ratio between low-depth kmers.
errorlowerconst=3   (elc) Use mult2 when the lower kmer is at most this deep.
mincountcorrect=3   (mcc) Don't correct to kmers with count under this.
pathsimilarityfraction=0.45(psf) Max difference ratio considered similar.
                           Controls whether a path appears to be continuous.
pathsimilarityconstant=3   (psc) Absolute differences below this are ignored.
errorextensionreassemble=5 (eer) Verify this many kmers before the error as
                           having similar depth, for reassemble.
errorextensionpincer=5     (eep) Verify this many additional bases after the
                           error as matching current bases, for pincer.
errorextensiontail=9       (eet) Verify additional bases before and after 
                           the error as matching current bases, for tail.
deadzone=0          (dz) Do not try to correct bases within this distance of
                    read ends.
window=12           (w) Length of window to use in reassemble mode.
windowcount=6       (wc) If more than this many errors are found within a
                    a window, halt correction in that direction.
qualsum=80          (qs) If the sum of the qualities of corrected bases within
                    a window exceeds this, halt correction in that direction.
rbi=t               (requirebidirectional) Require agreement from both 
                    directions when correcting errors in the middle part of 
                    the read using the reassemble algorithm.
errorpath=1         (ep) For debugging purposes.

Junk-removal parameters (to only remove junk, set mode=discard):
tossjunk=f          Remove reads that cannot be used for assembly.
                    This means they have no kmers above depth 1 (2 for paired
                    reads) and the outermost kmers cannot be extended.
                    Pairs are removed only if both reads fail.
tossdepth=-1        Remove reads containing kmers at or below this depth.
                    Pairs are removed if either read fails.
lowdepthfraction=0  (ldf) Require at least this fraction of kmers to be
                    low-depth to discard a read; range 0-1. 0 still
                    requires at least 1 low-depth kmer.
requirebothbad=f    (rbb) Only discard pairs if both reads are low-depth.
tossuncorrectable   (tu) Discard reads containing uncorrectable errors.
                    Requires error-correction to be enabled.

Shaving parameters:
shave=f             Remove dead ends (aka hair).
rinse=f             Remove bubbles.
wash=               Set shave and rinse at the same time.
maxshavedepth=1     (msd) Shave or rinse kmers at most this deep.
exploredist=300     (sed) Quit after exploring this far.
discardlength=150   (sdl) Discard shavings up to this long.
Note: Shave and rinse can produce substantially better assemblies
for low-depth data, but they are very slow for large metagenomes.
So they are recommended to be enabled for optimal results.

Overlap parameters (for overlapping paired-end reads only):
merge=f             Attempt to merge overlapping reads prior to 
                    kmer-counting, and again prior to correction.  Output
                    will still be unmerged pairs.
ecco=f              Error correct via overlap, but do not merge reads.
testmerge=t         Test kmer counts around the read merge junctions.  If
                    it appears that the merge created new errors, undo it.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.
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

z="-Xmx14g"
z2="-Xms14g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

tadpole() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP assemble.Tadpole $@"
	echo $CMD >&2
	eval $CMD
}

tadpole "$@"
