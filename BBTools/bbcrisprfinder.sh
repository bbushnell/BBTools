#!/bin/bash

usage(){
echo "
Written by Brian Bushnell with help from Simon Roux
Last modified Nov 22, 2023

Description:  Finds interspersed repeats contained within sequences;
specifically, only information within a sequence is used.  This is based
on the repeat-spacer model of crisprs.  Designed for short reads, 
but should work with full genomes.

Usage:  bbcrisprfinder.sh in=<input file>
See bbmap/pipelines/crisprPipeline.sh for more examples.

Standard parameters:
in=<file>       Primary input (a fasta or fastq file).
out=<file>      Output for annotated reads containing 
                interspersed repeats.  Reads not matching the criteria
                are discarded (though pairs are retained).
outc=<file>     Output for crispr sequences and their flanking repeats.
                Repeats are in lower-case; spacer is upper-case.
outr=<file>     Output of just repeats (one copy of each unique repeat).
                Only includes full-length repeats.
outs=<file>     Output of just spacers (one copy of each unique spacer).
                Only includes spacers adjacent to at least one full-length
                repeat (so the bounds are known).
chist=<file>    Histogram of crispr stats.
phist=<file>    Histogram of palindrome stats.
*All output files are optional.

Processing parameters:
merge=f         Try to merge paired reads before processing them.
masked=f        Indicates that the sequences are masked to lowercase;
                repeats will not be extended outside of lowercase regions.
                BBDuk or BBMask can produce masked sequences.
pad=0           Extend the boundaries of the repeats by this much when
                printing to outcrisper.  The padding will be in uppercase.
                This helps to show why a repeat was not extended.
annotate=f      Rename reads to indicate the repeat boundaries.  Annotations
                look like this: [56-77,117-138;150],P=7,M=5,L=6,S=62,T=6+0
                ...meaning the repeat occurs at positions 56-77 and 117-138;
                the read length is 150; and there is a palindrome of length
                7 with 5 matches, a loop length of 6, starting at position
                62, with tails outside the palindromic region length 6 and 0.
reads=-1        If positive, quit after processing this many reads.

Repeat Parameters:
krepeat=13      (kr) Use this seed kmer length to find repeats.  
                Repeats shorter than k may not be found.  Max is 31.
mm=1            (mmr) Mask this many bases in the middle of kmers to allow 
                imperfect repeats.  Also requires rmismatches to be set.
minrepeats=2    Ignore repeats with fewer than this many nearby copies.
                The more copies, the better consensus will refine the borders.
                For this purpose 2 partials count as 1 full-length repeat.
rmismatches=3   (rmm) Maximum allowed mismatches in a repeat pair.
minspacer=14    Ignore spacers shorter than this.
maxspacer=60    Ignore spacers longer than this.
minrepeat=20    Ignore repeats shorter than this.
maxrepeat=54    Ignore repeats longer than this.
minrgc=0.09     Ignore repeats with GC below this.
maxrgc=0.89     Ignore repeats with GC above this.
grow=t          Extend repeats through mismatches, if rmm>0.  Increases the
                number of unique repeats, and decreases unique spacers. 
lookahead=5     Minimum subsequent matches to extend through a mismatch.

Reference Parameters:
ref=<file>      If a reference of known CRISPR repeats is supplied, all
                detected repeates will be aligned to it, using the best match
                to override the predicted repeat boundaries.  Subsequent
                parameters in this section have no effect without a ref.
outref=<file>   Output the reference sequences used, annotated to indicate
                the number of uses and palindrome position.
kref=13         Kmer length for selecting ref sequences to attempt align.
                Lower is more sensitive but can take exponentially longer.
                'k' will set both krepeat and kref to the same value.
mmref=1         (maskMiddleRef) Number of bases masked in the middle of the
                kmer. 0 is umasked, and more increases mismatch tolerance.
minrefm=18      (minRefMatches) Reject alignments with fewer matching bases.
refmm=5         (refMismatches) Reject alignments with more mismatches.
refmmf=0.2      (refmismatchfraction) Allowed mismatches will be the maximum
                of refmm and refmmf*length.
minrefc=0       (minRefCount) Only load reference sequences with a count
                of at least this.  Intended for repeats generated via the
                'outr' flag, whose headers have a 'count=x' term indicating
                how many times that exact repeat was encountered.  Small
                numbers can be slow with large references.
discardaf=t     (discardAlignmentFailures) When the best alignment for a
                repeat fails repeat thresholds like minspacer, minrepeat, or
                rmismatches, discard that repeat.
shrinkaf=f      Trim mismatches from the ends of alignment failures.
revertaf=f      Revert alignment failures to their pre-alignment borders.
discardua=f     (discardUnaligned) Discard repeats that do not align to any
                ref repeats.  This means they did not meet minrefm/maxrefmm.
minrepeat0=11   Allow repeats this short prior to alignment (discarded unless
                they lengthen).  Has no impact if below kref or krepeat.
sortref=auto    Sort alignment candidates by count, then length.  May be set
                to t/f/auto. When set to auto, sequences will be sorted only
                if counts are present.
doublefetch=t   (ff) Fetch ref sequences using kmers from both repeat copies.
                Increases sensitivity in the presence of mismatches.
doublealign=t   (aa) Align ref sequences to both ref repeat copies.

Consensus Parameters:
consensus=t     When 3+ nearby repeat copies are found, adjust the boundaries
                so that all are identical.
minoverlapconsensus=18  Only try consensus on repeats that overlap at least
                        this much.
maxtrimconsensus=5      Do not trim more than this many bases on each end.

Partial Repeat Parameters:
bruteforce=t    Look for partial or inexact repeats adjacent to detected
                repeats.  These can improve consensus.
mintailrepeat=9 (mtr) Minimum length of partial repeats on read tips.
rmmt=1          Maximum allowed mismatches in short repeats at read tips.
rmmtpad=3       Ignore this many extra leading mismatches when looking
                for tip repeats (they will be trimmed later).

Palindrome Parameters:
minpal=5        If greater than 0, each repeat will be scanned for its 
                longest palindrome of at least this length (just the
                palindrome sequence, excluding the loop or tail).
                Palindromes will be annotated alongside repeats.
pmatches=4      Minimum number of matches in a palindrome.
pmismatches=2   Maximum allowed mismatches in a palindrome.
minloop=3       Ignore palindromes with a loop shorter than this.
maxloop=26      Ignore palindromes with a loop longer than this.
mintail=0       Ignore palindromes with a tail shorter than this.
maxtail=24      Ignore palindromes with a tail longer than this.
maxtaildif=21   Ignore palindromes with a tail length difference greater 
                than this.
reqpal=f        Discard repeats that lack a suitable palindrome.
symmetric=f     Trim repeats to make them symmetric around the palindrome.
                Not recommended.

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

crisprfinder() {
	local CMD="java $SIMD $EA $EOOM $z -cp $CP jgi.CrisprFinder $@"
	echo $CMD >&2
	eval $CMD
}

crisprfinder "$@"
