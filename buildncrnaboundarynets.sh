#!/bin/bash
# Learns and evaluates boundary corrections, then exports a loadable family.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 6, 2026

Usage: buildncrnaboundarynets.sh config=learning.config out=new_family_directory

Learns and compares classifier/regression boundary corrections, then retrains an
eligible winner on every eligible locus. If neither candidate, or its final
retraining, improves both ends without violating detection/precision limits,
the output family explicitly disables boundary refinement.
Use mode=refit to train an already chosen method once, without caller evaluation
or method selection. That export is labeled REFIT_UNEVALUATED.

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
familyroot=         Optional source family/group directory for manifest bindings.
                    Core ref/hbm/seeds still override its corresponding payloads.
                    Omission preserves installed-resource lookup. Export is checked
                    before training, including retained sibling/cutoff resources.
mode=build          Normal learn/evaluate/refit workflow.
mode=evaluate       Evaluate saved tuning products only; never train or refit.
learned=            Required with mode=evaluate: existing directory containing
                    classifier.bbnet, regression.bbnet, start.tsv and stop.tsv.
                    Use the original corpus, split, library and geometry settings.
                    Output is labeled TUNING_ONLY, not an all-data release refit.
mode=refit          Train one frozen method on all eligible loci and export it.
method=             Required only for refit: classifier or regression. Supply
                    only its classifierconfig= or regressionconfig= recipe.
                    Uses one training shard and monitors that same input.
                    No caller, selection or acceptance fallback runs; evaluate
                    the exported network with the intended inference policy.
                    Refit needs no truthfamily/genomes/truth/callerconfig,
                    splitseed or shards, and rejects those explicit controls,
                    plus learned/maxlost/minprecision/tiepreference and the
                    other trainer's config. No minimum TaxID-group count.
exclude=            Optional original locus IDs excluded only from learning.
shards=1            Length-balanced whole-locus training shards; one validation file.
features=v1         Ten shared inference features; currently the only feature set.
marginstart=0       Minimum score improvement needed to move the start.
marginstop=0        Minimum score improvement needed to move the stop.
cutoff=             Optional finite endpoint score cutoff; absent means no cutoff.
maxlost=0           Maximum baseline-detected reference losses.
minprecision=0.95   Minimum candidate precision, in [0,1]. Explicit 0 removes
                    the precision floor; maxlost and positive endpoint gain
                    still apply. The default remains 0.95.
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
Fixed refit records all source groups with withheld=0, table/vector counts and
any exclusions, missing context or model-selection failures. Its network is
loadable with refinement enabled; REFIT_UNEVALUATED is not release approval.

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
