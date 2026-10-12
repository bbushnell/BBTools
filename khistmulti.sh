#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Fischl
Last modified October 11, 2026

Description: Estimates kmer depth histograms for multiple kmer lengths in one
input pass. Each worker uses one compact DynamicDemiLog sketch per K; completed
sketches are merged. No full kmer tables or intermediate read files are needed.

Usage: khistmulti.sh in=reads.fq.gz out=hist.tsv k=32,64,96,124
       khistmulti.sh in1=R1.fq.gz in2=R2.fq.gz out=hist.tsv raw=t

Input/output:
in=<file>       FASTA/FASTQ/SAM/BAM input, interleaved pairs or single reads.
in2=<file>      Optional separate R2 file. Both mates contribute once.
int=auto        Detect interleaved FASTQ; int=t forces paired input.
out=stdout.txt  One tab-delimited histogram. khist= is an alias.
overwrite=t     Allow replacing an existing output. append is unsupported.
config=<file>   Read additional flag=value arguments from a text file.
reads=-1        Limit input fragments (reads or pairs); -1 means all.
qin=auto        Input quality encoding. minprob=0 ignores quality scores.

Histogram parameters:
k=32,64,96,124  Positive, unique K values in output-column order. Arbitrary K
                is supported, including K>read length (an empty histogram).
buckets=256k    Buckets per K per worker, rounded UP to a power of two.
                Larger sketches reduce sampling noise. Default actual=262144.
seed=0          Reproducible hash seed; a negative value chooses a random seed.
minprob=0       Exclude kmers below this probability of all bases being correct.
histmax=65535   Last bin includes ALL depths >=histmax; range 1..65535.
                DDL counters saturate at 65535 regardless of histmax.
raw=f          Add Raw_K columns containing sampled occupied-bucket counts.
zeros=f        Include rows whose counts are zero for every K.
t=auto         Worker threads. threadsin= controls input parsing threads.

Output: Depth rows and Count_K columns, estimated DISTINCT canonical kmers at
each observation depth, NOT read depth or genomic coverage. Counts are scaled
by cardinality/occupied buckets, with three decimal places. Raw counts are
sampled bucket winners, not an exact kmer census. DDL fingerprint ties can
combine counts. Sparse bins are noisy; compare distributions, not single bins.
The final bin is a lower bound. Metadata records seed, cardinalities and sample
sizes. No peak caller or automatic assembly-K selection is performed.

Memory: Sketch arrays use 5*buckets*number_of_Ks*workers bytes, plus reader,
object and histogram overhead. Default four Ks use 5 MiB per worker.

Java parameters:
-Xmx1g         Default maximum heap; increase for many Ks/buckets/threads.
-eoom          Exit on an out-of-memory error.
-da            Disable assertions (enabled by default).

Standard input: in=stdin.fq or in=stdin.fa. BGZF/gzip input is supported by
BBTools' native I/O. Histogram compression is selected by output extension.
Please contact Brian Bushnell at bbushnell@lbl.gov with problems.
"
}

if [ $# -eq 0 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit 0
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then CP="$DIR/bbtools.jar"; else CP="$DIR/current/"; fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"
	parseJavaArgs "--xmx=1g" "--xms=128m" "--mode=fixed" "$@"
	setEnvironment
}

resolveSymlinks
setEnv "$@"
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" cardinality.KmerHistMulti "$@"
