#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 29, 2026

Description:  Assigns taxonomy to query sequences by comparing kmer
frequencies to those in a reference database.  Developed for taxonomic
assignment of metagenomic bins, but it can also run on a per-sequence basis.
QuickClade is extremely fast and uses little memory.  However, the accuracy
declines for incomplete genomes.  k5dif represents the sum of the absolute 
values of the differences between the 5-mer frequency spectra, so 
the range is 0-1.  Because no marker genes are used, QuickClade should 
perform similarly for any clade in the reference dataset.
While the default reference is taxonomically labeled, you can use whatever
you want as a reference, with or without taxonomic labels.

Usage Examples:
quickclade.sh query1.fa query2.fa query3.fa
or
quickclade.sh bins
or
quickclade.sh contigs.fa percontig out=results.tsv usetree

For accuracy evaluation:
quickclade.sh printmetrics usetree genomesdir out=null includeself=f


File Parameters:
in=<file,file>  Query files or directories.  Loose file or directory names are
                also permitted.  Input can be fasta, fastq, or spectra files;
                spectra files are made by cladeloader.sh.
ref=<file,file> Reference files; the current default is:
                refseqA48_with_ribo.spectra.gz
                It is plaintext, human-readable, and pretty small.
out=stdout      Set to a file to redirect output.  Only the query results will
                be written here; progress messages will still go to stderr.
server          Use this flag to send kmer spectra to a remote server if you do not
                have a local database.

Presets (override individual settings; can be further overridden by later flags):
fast            records=1, buffer=1, callssu=f, sketch=f.  Fastest mode.
medium          records=5, buffer=20, callssu=t, sketch=t.  Default behavior.
slow            records=10, buffer=50, callssu=t, sketch=t, index=t.
                Automatically increases memory to 8g unless -Xmx is explicit.

Basic Parameters:
percontig       Run one query per contig instead of per file.
perfile         Opposite of percontig (default); one query per input file.
minlen=0        Ignore sequences shorter than this in percontig mode.
records=5       Print this many top hits per query.  Sets both cladehits
                and sketchhits.
cladehits=5     Max clade-based hits to display (independent of sketchhits).
sketchhits=5    Max sketch-based hits to display (only with sketch and index).
buffer=20       Internal candidate buffer size for ranking.  A larger buffer
                finds better top hits because SSU alignment is evaluated
                lazily.  Only affects clade hits, not sketch hits.
steps=6         Only search up to this many GC intervals (of 0.01) away from
                the query GC.
callssu=t       Call 16S and 18S for alignment to reference SSU.  Slightly
                slower.  Affects top hit ordering.
server=f        Send spectra to server instead of using a local reference.
                Enabled automatically if there is no local reference.
composition=    Output a taxonomy composition report to this file.  Shows
                per-level tables with bases, sequences, and percentages.
                Use composition=stdout to print to screen after results.
                Best used with percontig or multiple input files.
                summary= is an alias for composition=.
		
Output Format:
format=human    Output format.  Options: human (multi-line per hit),
                machine (one line per hit, tab-delimited header),
                tabular (compact one-line per hit with column headers).
showrecords=t   Set to false to suppress per-record output.  Useful with
                composition to display only the taxonomy summary.
color=t         ANSI color coding of hits by taxonomic level.  Default on
                for human and tabular formats, off for machine format.
colorlevel=     Taxonomic level for color grouping (default family).
                Hits in the same family get the same color; off-family
                hits are visually distinct.

showloading=t   Print loading progress messages to stderr.  Set to false
                to suppress index/sketch/query loading messages.

Taxonomy Filtering:
level=          Filter hits by taxonomic level (e.g. level=family).  Shows
                only the best hit per taxon at that level, so instead of
                7 E. coli strains you see the best hit from each family.
		Constrained to hits within the buffer.
topcount=10     Max entries per taxonomic level in the composition report.
minfraction=0   Minimum fraction (0-1) to include in composition report.

Proxy Parameters:
proxyhost=<addr>  HTTPS proxy hostname for environments requiring a proxy
                to reach external servers.  Sets -Dhttps.proxyHost for Java.
proxyport=<num>   HTTPS proxy port number.  Sets -Dhttps.proxyPort for Java.

DDL Sketch Parameters:
sketch=t        Enable sketch-based matching using DDL (DynamicDemiLog)
                cardinality profiles.  Loads refseqSketchDDL.tsv.gz from
                the resources directory by default.  Also builds a DDL
                from each query for comparison.  ddl=t is an alias.
sketchfile=     Path to a specific DDL sketch file.  Overrides the default.
                ddlfile= and sketchref= are aliases.
sketchindex=f   Build an index from DDL sketches; this allows hits by 31-mer
                matching, orthogonal to the clade index, allowing LCA
                (lowest common ancestor) calculation.  Implies sketch=t.
minsketchhits=3 Minimum matching DDL buckets to report a sketch hit.
ddlk=31         K-mer length for DDL sketches.
ddlbuckets=2048 Number of buckets in DDL sketches.

Threading Parameters:
loadthreads=auto  Number of threads for parsing reference spectra.
                  By default uses all available threads.
ddlloadthreads=auto  Number of threads for loading DDL sketch files.
                  DDL files load in parallel when multiple are present;
                  this controls the total thread budget across all files.
comparethreads=auto  Number of threads for query comparisons.
parallelsetup=t Load tree, reference index, and queries in parallel.
                Disable with parallelsetup=f for lower memory usage.

Advanced Parameters (mainly for benchmarking):
printmetrics    Output accuracy statistics; mainly useful for labeled data.
                Labeled data should have 'tid_1234' or similar in the header.
                Works best with 'usetree'.
printqtid       Print query TaxID.
banself         Ignore records with the same TaxID as the query.  Makes the
                program behave like that organism is not in the reference.
simd            Use vector instructions to accelerate comparisons.
maxk=5          Can be set to 4 or 3 to restrict kmer frequency comparisons
                to smaller kmers.  This may improve accuracy for small
                sequences/bins, but slightly reduces accuracy for large
                sequences/bins.
ccm=1.2         Threshold for using pentamers; lower is faster.
ccm2=1.6        Threshold for using tetramers.
gcdif=0.04      Initial maximum GC difference.
gcmult=0.5      Max GC difference as a fraction of best 5-mer difference.
strdif=0.12     Initial maximum strandedness difference.
strmult=1.2     Max strandedness difference as a fraction of best 5-mer diff.
hhdif=0.025     Maximum HH metric difference.
cagadif=0.017   Maximum CAGA metric difference.
hhmult=0.5      Max HH difference as a fraction of best 5-mer difference.
cagamult=0.8    Max CAGA difference as a fraction of best 5-mer difference.
ee=t            Early exit; increases speed.
entropy         Calculate entropy for queries.  Slow; negligible utility.
usetree         Load a taxonomic tree for better grading for labeled data.
aligner=quantum Options include ssa2, glocal, drifting, banded, crosscut.

Distance Metrics:
abs             Use absolute difference of kmer frequencies.
cos             Use 1-cosine similarity of kmer frequencies.
euc             Use Euclidian distance.
hel             Use Hellinger distance.
abscomp         GC-compensated version of abs (default).
Note:  The distance metric strongly impacts ccm, gcmult, and strmult.
       Defaults are optimized for abscomp.

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

	# Detect slow/index mode for RAM sizing
	local DEFAULTXMX="4g"
	local DEFAULTXMS="4g"
	local EXPLICIT_XMX=false
	for arg in "$@"; do
		case "$arg" in
			slow|sensitive) DEFAULTXMX="8g"; DEFAULTXMS="8g" ;;
			index|index=t|index=true|sketchindex|sketchindex=t|sketchindex=true)
				DEFAULTXMX="8g"; DEFAULTXMS="8g" ;;
			-Xmx*|--xmx=*) EXPLICIT_XMX=true ;;
		esac
	done
	if [ "$EXPLICIT_XMX" = true ]; then
		parseJavaArgs "--mode=fixed" "$@"
	else
		parseJavaArgs "--xmx=${DEFAULTXMX}" "--xms=${DEFAULTXMS}" "--mode=fixed" "$@"
	fi
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $PROXY $XMX $XMS -cp $CP clade.CladeSearcher $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
