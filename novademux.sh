#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 18, 2024

Description:  Demultiplexes sequencer reads into multiple files based on
their barcodes.  Uses statistical analysis to ensure optimal yield and
minimal crosstalk in the presence of errors.  Barcodes (indexes) must be
embedded in read headers, and the expected barcodes must be provided
as a text file with one barcode (or barcode pair) per line.

Usage:
novademux.sh in=reads.fq out=out_%.fq outu=unknown.fq expected=barcodes.txt

For twin files:
novademux.sh in=in_#.fq out=out_%_#.fq.gz outu=unk_#.fq expected=barcodes.txt


File Parameters:
in=<file>       Input file.
in2=<file>      If input reads are paired in twin files, use in2 for the 
                second file.  You can alternatively use the # symbol,
                e.g. 'in=read_#.fastq.gz', which is equivalent to
                'in1=read_1.fastq.gz in2=read_2.fastq.gz'.
out=<file>      Output files for reads with matched headers (must contain %
                symbol).  For example, out=out_%.fq with indexes XX and YY 
                would create out_XX.fq and out_YY.fq.  If twin files for
                paired reads are desired, use the # symbol.  For example,
                out=out_%_#.fq in this case would create out_XX_1.fq, 
                out_XX_2.fq, out_YY_1.fq, and out_YY_2.fq.
outu=<file>     Output file for reads with unmatched headers.
stats=<file>    Print statistics about how many reads went to each file.
expected=       List of barcodes (or files containing containing barcodes) 
                to parse from read headers.  Files should contain one barcode 
                per line.  For example, 'expected=barcodes.txt' or
                'expected=ACGTACGT,GGTTAACC,AACCGGTT'.  This list must be
                contain all pooled barcodes to ensure accuracy, including
                PhiX if present.
writeempty=t    Write empty files for expected but absent barcodes.
subset=         Optional list of barcodes when only some output files are
                desired; only demultiplex these libraries.
nosplit=f       When true, dump all reads to outu instead of individual files.
rename=f        When true, append the assigned barcode (or 'unknown') to
                each read header, after a tab.  Can be used in conjunction 
                with nosplit to simply label reads.
rc1=f           Reverse-complement index1 from expected and samplemap.
rc2=f           Reverse-complement index2 from expected and samplemap.
addpolyg=f      It is recommended to set this to true on a platform where
                no signal is read as G.  This will add poly-G as a dummy
                expected barcode.  If no signal yields a different base call,
                use the appropriate flag (addpolyc, etc).
remap=          Change symbols for output filenames.  For example, remap=+-
                would output barcode ACGT+TGCA to file ACGT-TCGA.fq.gz.

Legacy Output Stats File Support Parameters:
legacy=         Set this to a path like '.' to output legacy stats files.
samplemap=      An input csv or tsv containing barcodes and sample names,
                for legacy stats.  If present 'expected' can be omitted.
lane=0          Set this to a number to print the lane in legacy files.

Barcode Parsing Mode Parameters (choose one):
barcode         Parse the barcode automatically, assuming the standard
                Illumina header format.  This is the default.
header          Match the entire read header.
prefix          Match the prefix of the read header (length must be set).
suffix          Match the suffix of the read header (length must be set).
hdelimiter=     (headerdelimiter) Split the header using this delimiter,
                then select a term (column must be set).  Normally the 
                delimiter will be used as a literal string (a Java regular
                expression); for example, ':' or 'HISEQ'.  But there are
                some special delimiters which will be replaced by the symbol
                they name, because they can cause problems.
                These are provided for convenience due to OS conflicts:
                   space, tab, whitespace, pound, greaterthan, lessthan, 
                   equals, colon, semicolon, bang, and, quote, singlequote
                These are provided because they interfere with Java regular 
                expression syntax:
                   backslash, hat, dollar, dot, pipe, questionmark, star,
                   plus, openparen, closeparen, opensquare, opencurly
                In other words, to match '.', you should set 'hdelimiter=dot'.

length=0        For prefix or suffix mode, use this many characters from
                the read header.  Must be positive in these modes.
column=0        Select the term when using a header delimiter.  This is
                1-based (first term is column 1) so it must be positive.

Barcode Assignment Mode Parameters (choose one):
mode=prob       prob: Default mode.  Assigns reads to the bin where they most
                   likely belong, from gathering statistics across the pool.
                tile: Similar to prob, but calculates statistics on a per-tile
                   basis for higher precision.  This mode is recommended as
                   long as the tile numbers are in the read headers.
                hdist: Demultiplex reads to the bin with the fewest 
                   mismatches.  This is the fastest and least accurate mode.
                   Here, 'hdist' stands for 'Hamming distance'.
Note: prob and tile mode may require a license.

Server Parameters (for prob Mode only):
server=auto     true:  Barcode counts are sent to a remote server for 
                       processing, and barcode assignments are sent back.
                false: Barcode counts are processed locally.
                auto:  Sets flag to false unless the local machine contains
                       proprietary probabilistic processing code.

Sensitivity Cutoff Parameters for Prob/Tile Mode:
maxhdist=6     Maximum Hamming distance (number of mismatches) allowed.
                Lower values will reduce yield with little benefit.
pairhdist=f     When true, maxhdist will apply to the Hamming distance of
                both barcodes combined (if using dual indexes).  When false,
                maxhdist will apply to each barcode individually.
minratio=1m     Minimum ratio of probabilities allowed; k/m/b suffixes are
                allowed.  ratio=1m will only assign a barcode to a bin if
                it is at least 1 million times more likely to belong in that
                bin than in all others combined.  Lower values will increase
                yield but may increase crosstalk.
minprob=-5.6    Discard reads with a lower probability than this of belonging
                to any bin.  This is log10-scale, so -5 means 10^-5=0.00001.
                Lower values will increase yield but increase crosstalk.
                E.g., -6 would be lower than -5.
matrixthreads=1 More threads is faster but adds nondeterminism.
Note: These cutoffs are optimized for dual 10bp indexes.  For single 10bp
indexes, 'minratio=5000 minprob=-3.2' is recommended.

Sensitivity Cutoff Parameters for HDist Mode:
maxhdist=1      Maximum Hamming distance (number of mismatches) allowed.
                Lower values will reduce yield and decrease crosstalk.
                Setting maxhdist=0 will allow exact matches only.
pairhdist=f     When true, maxhdist will apply to the Hamming distance of
                both barcodes combined (if using dual indexes).  When false,
                maxhdist will apply to each barcode individually.
clearzone=1     (cz) Minimum difference between the closest and second-closest
                Hamming distances.  For example, AAAG is 1 mismatch from
                AAAA and 3 mismatches away from GGGG, for a margin of 3-1=2.
                This would be demultiplexed into AAAA as long as the
                clearzone is set to at most 2.  Lower values increase both
                yield and crosstalk.

Buffering Parameters:
streams=8       Allow at most this many active streams.  The actual number
                of open files will be 1 greater than this if outu is set,
                and doubled if output is paired and written in twin files 
                instead of interleaved.  Setting this to at least the number
                of expected output files can make things go much faster.
minreads=0      Don't create a file for fewer than this many reads; instead,
                send them to unknown.  This option will incur additional
                memory usage.
rpb=8000        Dump buffers to files when they fill with this many reads.
                Higher can be faster; lower uses less memory.
bpb=8000000     Dump buffers to files when they contain this many bytes.
                Higher can be faster; lower uses less memory.

Spike-in Processing Parameters (particularly for spike-ins with no barcodes):
spikelabel=     If and only if a spike-in label is set here, reads will be
                aligned to a reference, and matching reads will be sent to
                the file with this label.  May be a barcode or other string.
refpath=phix    Override this with a file path for a custom reference.
kspike=27       Use this kmer length to map reads.
minid=0.7       Identity cutoff for matching the reference.
mapall=f        Map all reads to the reference, instead of just unassigned
                reads.

Common Parameters:
ow=t            (overwrite) Overwrites files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=auto        (interleaved) Determines whether INPUT file is considered 
                interleaved.                

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify
                200 megs.  The max is typically 85% of physical memory.
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

z="-Xmx16g"
z2="-Xms16g"
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
	freeRam 16000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function demux() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP barcode.NovaDemux $@"
	echo $CMD >&2
	eval $CMD
}

demux "$@"
