#!/bin/bash

usage(){
echo "
Written by Noire and Brian Bushnell
Last modified August 24, 2026

Description:  Calibration and coverage measurement tool for ncRNA consensus
libraries.  A mini CallGenes for calibrating each individual ncRNA class.

For each query, shortlists reference models via a 7-mer inverted index
(same infrastructure as the scavenger), aligns in decreasing kmer-sharing
order, and reports mapping/alignment pass rates plus diagnostic histograms.

Output follows BBDuk conventions: out= captures queries that FAIL to match,
outm= captures queries that PASS (mapped and above identity threshold).

Usage:  ncrnacalibrate.sh in=<queries.fa> ref=<consensus.fa> out=<uncovered.fa>

File parameters:
in=<file>       Query sequences (e.g. genes extracted via cutgff with flank=20).
ref=<file>      Reference consensus library (e.g. rnasep_consensus.fa).
out=<file>      Output: queries that FAIL (unmapped or below minid).
outm=<file>     Output: queries that PASS (mapped and above minid).

Calibration parameters:
trim=0          Trim this many bases from each end before output
                (e.g. trim=20 to remove flanking added by cutgff flank=20).
minid=0.75      Minimum alignment identity to pass.
minkmers=0      Minimum shared 7-mers with top model to pass mapping.
aligner=quantum Choice of aligner: quantum, scrabble, or glocal.

Index parameters (same as scavenger, tune per family):
indexk=7        Kmer length for the model-selection index.
topn=60         Maximum models on the shortlist per query.
adaptive=t      Use adaptive shortlist cutoff (vs fixed).
floor=11        Adaptive: minimum shared kmers to make shortlist.
topfrac=0.48    Adaptive: minHits = this fraction of best model's count.
qfrac=0.072     Adaptive: minHits = this fraction of query kmer count.
fixedminhits=12 Fixed mode: minimum shared kmers.

Histogram output (default stderr; set to a filename to write to file):
anihist=stderr  Best-identity-per-query histogram.
kmerhist=stderr Shared-kmer-count histograms (top hit and accepted hit).
ratiohist=stderr  Kmer ratio (accepted/top) histogram.
histbins=100    Number of bins for the ANI histogram.

Output reports:
  Passed mapping:   N  X.XXX%   (of total queries)
  Failed mapping:   N  X.XXX%   (of total queries)
  Passed alignment: N  X.XXX%   (of passed mapping)
  Failed alignment: N  X.XXX%   (of passed mapping)

Java Parameters:
-Xmx            Set memory usage.  Default autodetected.

Please contact Brian at bbushnell@lbl.gov if you encounter any problems.
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

resolveSymlinks(){
	if [ -d "$DIR""current/" ]; then
		CP="$DIR""current/"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.NcrnaCalibrator $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
