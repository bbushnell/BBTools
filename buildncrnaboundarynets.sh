#!/bin/bash
# Learns and evaluates boundary corrections, then exports a loadable family.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 1, 2026

Usage: buildncrnaboundarynets.sh config=learning.config out=new_family_directory

Learns and compares classifier/regression boundary corrections, then retrains an
eligible winner on every eligible locus. If neither candidate, or its final
retraining, improves both ends without violating detection/precision limits,
the output family explicitly disables boundary refinement.

Required learning inputs:
in=                 Sense-oriented flanked FASTA; exact locus IDs with source=,
                    lflank=, rflank= and whole-contig contig_gc= metadata.
roster=             Complete original genome roster TSV.
loci=               Complete original annotation ledger TSV.
graph=              Version-two ORIGINAL original.support.tsv from
                    reducencrnamodels.sh; the old buildsupportgraph.sh output
                    uses a different format and is not accepted here.
original=           Original consensus FASTA matching the graph model axis.
ref=                Accepted consensus FASTA for vector model choice and calling.
hbm=                Accepted index-aligned native HBM/HBMC library.
seeds=              Accepted forward 17-mer FASTA.
family=             One ordinary shared-boundary family member, e.g. sam1.
truthfamily=        Exact label in the evaluation truth table.
genomes=            Text file of distinct absolute evaluation genome FASTA paths.
truth=              TSV: seqid, family, start, stop, strand; one-based inclusive.
callerconfig=       Ordinary combined CallGenes policy; caller defaults still apply.
                    No input/output, selected resource/refinement or selected-family
                    mapping overrides. Registered mapping settings train and call.
classifierconfig=   Native train.sh recipe, dims=10,...,1, explicit seed and duration,
                    vfraction=0, one network/cycle. No data/output routing flags.
regressionconfig=   Native regressiontrainer.sh recipe with the same constraints.
splitseed=          Nonnegative whole-TaxID split seed.
startradius=        Nonnegative search radius, no greater than library median length.
stopradius=         Same for the 3-prime end; at least one end needs negative examples.
startinside=        Nonnegative bases inside the start endpoint word.
startoutside=       Nonnegative bases outside; inside+outside must be 1..15.
stopinside=         Corresponding stop-end geometry.
stopoutside=        Corresponding stop-end geometry.
out=                Fresh final family directory; its parent must exist.

Options:
exclude=            Optional original locus IDs excluded only from learning.
shards=1            Length-balanced whole-locus training shards; one validation file.
features=v1         Ten shared inference features; currently the only feature set.
marginstart=0       Minimum score improvement needed to move the start.
marginstop=0        Minimum score improvement needed to move the stop.
cutoff=             Optional finite endpoint score cutoff; absent means no cutoff.
maxlost=0           Maximum baseline-detected reference losses.
minprecision=0.95   Minimum candidate precision, in [0.95,1].
tiepreference=regression  Final tie choice: classifier or regression.
t=auto             Worker budget; up to three caller trials run concurrently.
-Xmx4g             Parent heap. Aggregate concurrent child heaps use the same
                    budget IN ADDITION to the parent. Request memory accordingly.

out.work must also be fresh. It retains split membership, tuning/refit tables,
vectors, networks, ordinary caller GFFs/configurations and comparison reports.
out itself is loadable with callgenes.sh ncrnafamily=... ncrnafamilyroot=out.
Grouped families export the whole group while preserving siblings and shared seeds.
No resources are installed or staged. Final retraining consumes former validation
loci: its evaluation is coverage evidence, not independent holdout performance.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
while [ -h "$SCRIPT" ]; do
	DIR="$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$SCRIPT")"
	SCRIPT="$RESOLVED_SCRIPT"
	if [ "${SCRIPT#/}" = "$SCRIPT" ]; then
		RESOLVED_SCRIPT="$DIR/$SCRIPT"
		SCRIPT="$RESOLVED_SCRIPT"
	fi
done
DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.BuildNcrnaBoundaryNets "$@"
