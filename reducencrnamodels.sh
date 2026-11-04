#!/bin/bash
# Proposes a smaller consensus library with optional caller protection and regeneration.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 1, 2026

Usage: reducencrnamodels.sh in=genes.fa ref=consensus.fa cutoff=0.70 out=new_directory

Proposes a greedy subset of ncRNA models at a declared identity threshold.
Input genes must be complete, sense-oriented FASTA records with distinct IDs.
Quantum aligns every gene to every model. The original and reduced support
graphs retain every gene, including genes with no supporting model.

Required:
in=             Gene FASTA, optionally compressed.
ref=            Consensus FASTA.
cutoff=         Inclusive identity threshold in [0,1]; no default.
out=            New directory beneath an existing parent; existing output refused.

Options:
hbm=            Index-aligned HBM/HBMC input; subset preserves its codec.
target=0.98     Fraction of ALL genes that must retain a supporter, in (0,1].
t=auto          Alignment workers, using the standard BBTools thread option.

Optional caller protection (supply all together, including hbm):
genomes=        Text file listing distinct absolute genome FASTA paths.
truth=          TSV: seqid, family, start, stop, strand; one-based inclusive.
family=         CallGenes family key, such as sam1.
truthfamily=    Exact truth-table label, such as RF00162.
callerconfig=   Ordinary CallGenes scientific options, including family activation
                and unchanged seed resources; no input/output or JVM options.
callerheap=2g   Heap for one sequential child caller, in addition to this JVM.
maxcalls=4      Maximum full caller passes, including the original baseline; at least 2.

Optional serial regeneration (requires caller protection):
regenerate=f   Rebuild selected models from their current members, in order.
consensusconfig=  Public consensus recipe: passes, minid, endtrim, stochastic,
                stochasticseed. The stage owns grouping, files and threads.
consensusheap=1g  Heap for the single consensus child, in addition to this JVM.
maxregencalls=5  Maximum NEW caller passes for regeneration, separate from
                maxcalls; at least 1. The verified reduction result is reused as baseline.

Regeneration preserves model identity/count and the genes supported on entry to
that stage. Greedy target coverage remains a separate measurement against the
original corpus. Failed member proposals leave both model and assignments intact.
Lost caller detections restore the baseline winner's complete consensus/HBM pair.
On no progress, budget exhaustion or final corpus-support loss, export the
verified pre-regeneration library and report that fallback explicitly.

Protection preserves baseline-detected reference loci, not raw baseline calls.
It currently supports families with CallGenes' dynamic paired-model overrides;
fixed rRNA and tRNA model paths are not supported by this protection interface.
Lost loci restore their original winning models as a batch, followed by a new
caller pass. If restoration stalls or the call budget expires, retain the
verified original library. Without genomes/truth, output is only a proposal.
No global-minimum or false-positive-rate claim is made.

Outputs: consensus.fa, optional models.hbm/hbmc, original.support.tsv,
members.tsv (reduced graph), summary.txt. Protected runs also retain full GFFs,
executed caller configurations, caller_results.tsv and lost-reference tables.
Regeneration adds the regeneration/ directory with models.tsv,
series.assignments.tsv, candidate.assignments.tsv, final.assignments.tsv,
caller_results.tsv (including strand-aware exact ends), caller_restorations.tsv,
exported.caller.tsv, exported.models.tsv, summary.tsv and native child outputs. Rejected proposals remain diagnostic, never the top-level library.
The output is a library subset, not an installable family resource directory.

Standard BBTools options use flag=value; config=, -Xmx and outstream= work normally.

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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.ReduceNcrnaModels "$@"
