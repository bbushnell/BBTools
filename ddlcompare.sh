#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Ady
Last modified August 3, 2026

Description:  Pairwise genome comparison using DynamicDemiLog (DDL) bucket
matching.  Creates a DDL sketch for each input, compares them, and reports
WKID, ANI, cardinality, containment, completeness, and bucket-level
statistics.  Can also compare a query against a pre-built DDL reference file,
or run collision tests on a DDL file.

Usage:  ddlcompare.sh genome1.fa genome2.fa
    or: ddlcompare.sh query.fa ref=ddls.tsv records=10
    or: ddlcompare.sh query.fa refseq records=10
    or: ddlcompare.sh qf=queries.tsv ref=ddls.tsv t=32
    or: ddlcompare.sh ref=ddls.tsv collisiontest

Pairwise Parameters:
in=<file>       First input file (query).  Also accepts positional arguments.
in2=<file>      Second input file (reference).

Reference Mode Parameters:
ref=<file>      Pre-built DDL reference file (TSV format from DDLLoader).
refseq          Shorthand for ref=resources/refseqSketchDDL_k25e5b4096.tsv.gz.
                Also accepts ref=refseq.
queryfile=<file> Pre-built DDL query file (TSV format).  Compares all queries
qf=<file>       against all references (multi-query batch mode).
records=20      Max hits to display.
minhits=5       Minimum matching DDL buckets to report a hit.
index=f         Use inverted index for query acceleration.
csr=t           Inverted-index storage: CSR packed arrays (default, far less RAM)
                or the matrix reference (csr=f).  Results are identical either way.
t=auto          Number of threads (default: all available cores).

Collision Test:
collisiontest   Measure all-pairs collision rate in a DDL reference file.
                Requires ref= to be set.

Blacklist:
blacklist=<file>  Kmer blacklist to apply while sketching queries.  If unset,
                autoloads refseqGenomeDDLBlacklist_k25e5b65536_fused.fa.gz from
                resources/ if present.  It is a kmer set, so it masks any bucket
                count.  Filters taxonomically uninformative kmers during query
                sketch construction; pre-built reference sketches already have
                blacklisting baked in.

Other Parameters:
ssu=f           Load SSU sequences and align them for the reported hits.
percontig=f     One sketch per input sequence rather than per file.
minsketch=400   Skip sequences shorter than this in percontig mode.
pload=t         Load SSU maps in parallel with the reference sketches.

Sketch Parameters:
k=25            K-mer length for hashing.
buckets=2048    Number of DDL buckets.
exponent=5      Exponent bits (1-8).  A reference or query file's #exponent
                header overrides this; files that disagree are a fatal error.

Examples:
ddlcompare.sh ecoli.fa mruber.fa
ddlcompare.sh ref.fa mutant.fa k=25 buckets=2048
ddlcompare.sh query.fq.gz ref=refseqSketchDDL_k25e5b2048.tsv.gz records=5
ddlcompare.sh qf=top1k.tsv.gz ref=refseqSketchDDL_k25e5b2048.tsv.gz t=32

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLCompare $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
