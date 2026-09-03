#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 21, 2026

Description:  Error corrects reads and/or filters by depth, storing
kmer counts in a count-min sketch (a Bloom filter variant).
This uses a fixed amount of memory.  The error-correction algorithm is taken
from Tadpole; with plenty of memory, the behavior is almost identical to 
Tadpole.  As the number of unique kmers in a dataset increases, the accuracy 
decreases such that it will make fewer corrections.  It is still capable
of making useful corrections far past the point where Tadpole would crash
by running out of memory, even with the prefilter flag.  But if there is
sufficient memory to use Tadpole, then Tadpole is more desirable.

Because accuracy declines with an increasing number of unique kmers, it can
be useful with very large datasets to run this in 2 passes, with the first
pass for filtering only using a 2-bit filter with the flags tossjunk=t and
ecc=f (and possibly mincount=2 and hcf=0.4), and the second pass using a
4-bit filter for the actual error correction.  Set passes=2 to do this
internally in a single invocation; see the Multipass parameters below.

Usage:  bbcms.sh in=<input file> out=<output> outb=<reads failing filters>

Example of use in error correction:
bbcms.sh in=reads.fq out=ecc.fq bits=4 hashes=3 k=31 merge

Example of use in depth filtering:
bbcms.sh in=reads.fq out=high.fq outb=low.fq k=31 mincount=2 ecc=f hcf=0.4

Error correction and depth filtering can be done simultaneously.

Multipass parameters:
passes=1        Set to 2 to run an internal 2-pass pipeline: pass 1 builds a
                cheap K=31 2-bit filter and discards junk reads (tossjunk=t,
                mincount=2, hcf=0.4, ecc=f), then pass 2 runs your full
                command line (k=, bits=, ecc=, etc., unchanged) on the
                survivors.  Pass 2 is the normal pass; passes= only controls
                whether a cheap pre-filtering pass runs first.  Recommended
                for huge, low-depth metagenomes with heavy depth-1 junk,
                where it is both faster and more accurate than one pass.

File parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary read output.
out2=<file>     Read 2 output if reads are in two files.
outb=<file>     (outbad/outlow) Output for reads failing mincount.
outb2=<file>    (outbad2/outlow2) Read 2 output if reads are in two files.
extra=<file>    Additional comma-delimited files for generating kmer counts.
ref=<file>      If ref is set, then only files in the ref list will be used
                for kmer counts, and the input files will NOT be used for
                counts; they will just be filtered or corrected.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Hashing parameters:
k=31            Kmer length.  Any value 1 and up is supported; k>31 uses a
                slower multi-word kmer representation, exactly representable
                by default (see packed= below); set packed=f to fall back
                to the legacy even-split layout, which silently rounds some
                k>31 values down to the nearest achievable length (a note is
                printed to stderr when this happens).
packed=t        For k>31, use a fully-packed multi-word kmer representation
                (non-symmetric: leading words are completely full and only
                the last word is partial).  Allows any k>31 to be
                represented exactly, with no rounding.  Slightly slower.
                Set to f for the legacy even-split layout.
fullmix=t       For k>31, fully mix every word of the kmer, including the
                first, into the hash function; improves Bloom-filter
                collision resistance at k>31.  No effect at k<=31.  Set to
                f to disable (not recommended in combination with packed=t).
hashes=3        Number of hashes per kmer.  Higher generally reduces
                false positives at the expense of speed; rapidly
                diminishing returns above 4.
dualhash=t      Use two independent hashes to select Bloom lanes.  Every lane,
                including the first, still passes through the normal mixer.
                For k>31, uses xor1 plus xor2; for shorter kmers, derives an
                independent secondary hash from the exact packed key.
                Enabled by default because BBCMS targets metagenomic datasets
                large enough for primary 64-bit hash collisions to matter.
                Set to f to recover the legacy single-key lane recurrence.
ksmall=         Optional sub-kmer length; setting to slightly lower than k
                can improve memory efficiency by reducing the number of hashes
                needed.  e.g. 'k=31 ksmall=29 hashes=2' has better speed and
                accuracy than 'k=31 hashes=3' when the filter is very full.
                Not supported when k>31; ignored (forced equal to k) in
                that case.
minprob=0.5     Ignore kmers with probability of being correct below this.
memmult=1.0     Fraction of free memory to use for Bloom filter.  1.0 should
                generally work; if the program crashes with an out of memory
                error, set this lower.  You may be able to increase accuracy
                by setting it slightly higher.
cells=          Option to set the number of cells manually.  By default this
                will be autoset to use all available memory.  The only reason
                to set this is to ensure deterministic output.
seed=0          This will change the hash function used.  Useful if running
                iteratively with a very full table.  -1 uses a random seed.
symmetricwrite=t  (sw) Increases counting accuracy for a slight speed penalty.
                Could be slow on very low-complexity sequence.
                
Depth filtering parameters:
mincount=0      If positive, reads with kmer counts below mincount will
                be discarded (sent to outb).
hcf=1.0         (highcountfraction) Fraction of kmers that must be at least
                mincount to pass.
requireboth=t   Require both reads in a pair to pass in order to go to out.
                When true, if either read has a count below mincount, both
                reads in the pair will go to outb.  When false, reads will
                only go to outb if both fail.
tossjunk=f      Remove reads or pairs with outermost kmer depth below 2.
(Suggested params for huge metagenomes: mincount=2 hcf=0.4 tossjunk=t)

Error correction parameters:
ecc=t           Perform error correction.
markerrors=f    Mark bounded low-count kmer runs as N instead of correcting.
                Requires two high-depth kmers on both sides; exact K and K-1
                runs are marked, while end/merged runs are reported and skipped.
bits=           Bits used to store kmer counts; max count is 2^bits-1.
                Supports 2, 4, 8, 16, or 32.  16 is best for high-depth data;
                2 or 4 are for huge, low-depth metagenomes that saturate the 
                bloom filter otherwise.  Generally 4 bits is recommended for 
                error-correction and 2 bits is recommended for filtering only.
ecco=f          Error-correct paired reads by overlap prior to kmer-counting.
merge=t         Merge paired reads by overlap prior to kmer-counting, and 
                again prior to correction.  Output will still be unmerged.
smooth=3        Remove spikes from kmer counts due to hash collisions.
                The number is the max width of peaks to be smoothed; range is
                0-3 (3 is most aggressive; 0 disables smoothing).
                This also affects tossjunk.
                

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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP bloom.BloomFilterCorrectorMultipass $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
