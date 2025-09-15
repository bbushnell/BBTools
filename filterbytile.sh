#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 18, 2024

Description:  Filters reads based on positional quality over a flowcell.
Quality is estimated based on quality scores and kmer uniqueness; ideally,
the input to this program will already be adapter-trimmed and 
quality-score recalibrated.  PhiX, if present, will be used to estimate
absolute error rates from kmer uniqueness rates.
All reads within a small unit of area called a micro-tile are averaged,
then the micro-tile is either retained or discarded as a unit depending
on whether its metrics fall outside of a specified range (generally some
number of standard deviations away from average).
This program can process individual libraries, but achieves optimal 
performance when processing one lane at a time.  For this purpose it is best
to use as much memory as possible (e.g., 200GB RAM for 5 billion reads),
though it will still work with much less memory.

Please read bbmap/docs/guides/FilterByTileGuide.txt for more information.


Usage:	filterbytile.sh in=<input> out=<output>

Input parameters:
in=<file>           Primary input file.
in2=<file>          Second input file for paired reads in twin files.
indump=<file>       Specify an already-made dump file to use instead of
                    analyzing the input reads.
barcodes=<file>     Optional list of expected barcodes, one per line.
reads=-1            Process this number of reads, then quit (-1 means all).
interleaved=auto    Set true/false to override autodetection of the
                    input file as paired interleaved.

Output parameters:
out=<file>          Output file for filtered reads.
dump=<file>         Write a summary of quality information by coordinates.
                    This can be later used for filtering individual libraries.
counts=<file>       Write barcode counts.


Tile parameters:
xsize=500           Initial width of micro-tiles.  For NovaSeqX use 520.
ysize=500           Initial height of micro-tiles.  For NovaSeqX use 590.
size=               Allows setting xsize and ysize to the same value.
target=1600         Iteratively widen the micro-tiles until they contain
                    an average of at least this many reads.
alignedreads=250    Average aligned reads per tile for error rate calibration.

A micro-tile is discarded if any of several metrics indicate a problem.
The metrics are kmer uniqueness (u), average quality (q), probability
of being error-free (e), and poly-G rate (pg).  
Each has 3 parameters: deviations (d), fraction (f), and absolute (a).  
After calculating the difference (delta) between a micro-tile and average, 
it is discarded only if all three of these conditions are true for at least
one metric (using quality as the example):
1) delta is greater than (qd) standard deviations.
2) delta is greater than average times the fraction (qf).
3) delta is greater than the absolute value (qa).
Tiles are also marked for discard if they have too few reads to calculate
statistics or an inferred error rate (ier) above an absolute value; ier
does not need deviations because it is calibrated. 

Filtering parameters:
udeviations=1.5     (ud) Standard deviations for uniqueness discarding.
qdeviations=2.4     (qd) Standard deviations for quality discarding.
edeviations=3.0     (ed) Standard deviations for error-free probablity. 
pgdeviations=1.4    (pgd) Standard deviations for poly-G discarding.
ufraction=0.01      (uf) Min fraction for uniqueness discarding.
qfraction=0.08      (qf) Min fraction for quality discarding.
efraction=0.2       (ef) Min fraction for error-free probablity discarding.
pgfraction=0.2      (pgf) Min fraction for poly-G discarding.
uabsolute=1         (ua) Min absolute value for uniqueness discarding.
qabsolute=2.0       (qa) Min absolute value for quality discarding.
eabsolute=6         (ea) Min absolute value for error-free probablity.
pgabsolute=0.2      (pga) Min absolute value for poly-G discarding.
ier=0.012           (inferrederrorrate) Maximum predicted base error rate.
                    A more recent addition and usually superior to using
                    uniqueness deviations, if ~1% PhiX is spiked in.
mdf=0.4             (maxdiscardfraction) Don't discard more than this 
                    fraction of tiles no matter how bad the data is.

Alignment parameters:
Note: Alignment will only be performed if there is no input sam file,
and nothing will go to the output sam file unless internal alignment occurs.
samin=<file>        Optional aligned sam input file for error rate analysis.
samout=<file>       Output file for aligned reads.  Can be sam or fastq.
align=true          If no sam file is present, align reads to the reference.
alignref=phix       Reference for aligning reads if there is no sam file.
alignk1=17          Kmer length for seeding alignment to reference.
alignk2=13          Kmer length for seeding alignment of unaligned reads
                    with an aligned mate.
minid1=0.62         Minimum identity to accept individual alignments.
minid2=0.54         Minimum identity for aligning reads with aligned mates.
alignmm1=1          Middle mask length for alignk1.
alignmm2=1          Middle mask length for alignk2.

Note: Alignment is optional, but allows translation of kmer depth to error 
rate at high resolution.  The default reference, phiX, is nonrepetitive down
to k=13.  For internal alignment, the reference must be a short single 
sequence that is almost completely nonrepetitive at the selected kmer length.
If a sam file is used, any reference is OK and the alignment parameters are
ignored, but it should have few mutations.  A SNP rate of 1/1000 (like human)
is acceptable but sets an inferred error rate floor of 0.001 (Q30).

Other parameters:
usekmers=t          Load kmers to calculate uniqueness and depth.
lowqualityonly=t    (lqo) Only discard low quality reads within bad areas, 
                    rather than the whole micro-tile.  This usually discards
                    most of the reads in the bad micro-tiles anyway.
recalibrate=f       Recalibrate reads while filling tile info.
                    Requires calibration matrices from CalcTrueQuality.
                    Changes sam output, but not positionally-filter output.
dmult=-.1           Lower increases amount removed when lqo=t.  At 0, only 
                    reads with below average quality (or polyG) are removed.
idmaskwrite=15      A bitmask, (2^N-1), controlling the fraction of kmers
                    loaded in the bloom filter.  15 means 1/16th are loaded.
                    0 uses all kmers.
idmaskread=7        Controls fraction of kmers read when counting uniqueness.
k=31                Kmer length for Bloom filter (uniqueness calculation).
hashes=3            Bloom filter hashes.
cbits=2             Bloom filter bits per cell.
merge=f             Merge reads for insert and error rate statistics.
                    This can make the program take ~50% longer and only
                    affects the dump file.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 GB of RAM; -Xmx200m will specify 
                    200 MB.  The max is typically 85% of physical memory.
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

z="-Xmx8g"
z2="-Xms8g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

filterbytile() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP hiseq.AnalyzeFlowCell $@"
	echo $CMD >&2
	eval $CMD
}

filterbytile "$@"
