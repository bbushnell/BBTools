#!/usr/bin/env bash
# Member-pooled positional tables and native train/monitor matrices; no training.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Build member-pooled single-table positional inputs from every eligible original locus.
Usage: buildncrnapositionalinputs.sh in=contexts.fa out=fresh_directory \
 family=registered_name familyroot=working_family ref=consensus.fa hbm=models.hbmc \
 seeds=seeds.fa.gz original=originals.fa graph=support.tsv roster=genomes.tsv \
 loci=loci.tsv histogram=member.histogram.tsv radii=6,10,16 samples=5 \
 trainseed=74003 monitorseed=74004 prior=1 t=16

All scientific arguments above are required. exclude= optionally names an explicit
locus exclusion file. No random subsampling, limits, overwrite or append.
teachingroutes= and teachingmodels= optionally restrict full-error split-trained
teaching to frozen TIP model routes using a complete saved vectors.loci.tsv.
This restriction is separate from outer exclusions. Fresh training winners must
match the saved native Quantum model axis; ordinary inputs are unchanged without it.
Every eligible geometry-complete locus teaches each table. Monitoring resamples
those same eligible loci and is NOT independent validation. Signed jitter draws
only observed errors within each radius; excluded tails are reported. Empty windows
fail with accounting and no invented fallback. prior= is TOTAL mass per position.
startinside=9 startoutside=0 stopinside=9 stopoutside=0
                 Semantic placement per end; each word length must be5..9.
                 START begins at endpoint-startoutside; STOP begins at
                 endpoint-stopinside+1. Defaults retain fully-inside K9.
                 Per-key smoothing is prior/4^K for each end's actual K.
                 Larger word lengths require a separate storage design.
tripletables=f   Optional equal-K7..9 experiment: add endpoint-exclusive left(-K) and
                 right(+1) tables beside the declared placement. Uses leave-one-
                 locus-out table scores for BOTH training and same-locus monitor.
                 Adds sliding-defined selected-HBM endpoint coverage at -1/0/+1
                 to the usual three extras. Width is 3*(2R+1)+6; at R10,69 inputs
                 and21 output classes. Fuzziness uses the jittered/unmoved span,
                 never truth coordinates; identity remains Quantum. Default path
                 and its ordinary single-table output bytes are unchanged.
interchange=f    Optional plain TSV delivery for tripletables=t radii=10.
evaluationseed=74005
                 Required third distinct draw seed with interchange=t. Explicit
                 exclude= loci form a separate held-out evaluation set; partial
                 loci stay excluded. Training/monitoring use LOO; evaluation
                 teaches no table and subtracts nothing. Any exact assembly
                 shared across training/evaluation fails. Missing TaxIDs are OK.
                 Six matrices each start with #dims TAB69 TAB21 TAB0, followed
                 only by90 numeric columns. No second matrix copy is written.
                 Matching .rows.tsv sidecars retain locus/source/assembly/split,
                 model/sample and zero-based inclusive raw/true sense coordinates.
fullerrors=f    With interchange=t, use the COMPLETE signed histogram, including
                 its measured zero and distant-error mass. No clamping/redraw.
                 Each probe uses the observed NONZERO signed-error weights on
                 separately normalized offset frequencies, with exact locus LOO
                 for training/monitoring and no subtraction for evaluation.
                 Radius remains10; headers become #dims TAB69 TAB22 TAB0 and
                 each row has91 numbers. Target22 is out-of-range truth (index21).
                 Writes explicitly labeled weighted-score tables and directional
                 error bands. Invalid spans and clipped/ambiguous words remain
                 counted limitations, never out-of-range training rows. These
                 experimental tables are not legacy caller resources.
enrichmentbits=0 With fullerrors=t, enrichmentbits=12 applies unsigned 6+6
                 ratio quantization after exact whole-locus subtraction for
                 training/monitoring. Evaluation/serving use cached decoded
                 float-log values. All three probes keep their own weighted
                 backgrounds. Writes dense native A48 tables as weighted.tsv.gz,
                 plus tables.bytes.tsv before vector generation. Omission keeps
                 existing score/vector bytes; no evidence pruning is performed.
refit=f         Explicit all-data product mode, requiring fullerrors=t.
                 Omit exclude=, evaluationseed= and fixed evaluation-row/offset options.
                 Uses every complete eligible original locus, including former
                 held-outs. Writes only train and same-locus monitor matrices
                 with the same LOO 69-input/22-output schema; no evaluation set
                 or independence claim. Learned residual histograms must also
                 come from the all-data population.
evaluationrows5=existing.5prime.rows.tsv
evaluationrows3=existing.3prime.rows.tsv
modeloffsets=training_supported.offsets.tsv
offsetminimum=10
evaluationcounts=accepted_start_count,accepted_stop_count
                 Optional fixed-trial centering with fullerrors=t; supply all
                 five bindings together. The histogram must be the separately
                 prepared training-only residual distribution. All six weighted
                 tables are rebuilt; training/monitor draws remain unchanged.
                 Evaluation reuses the exact original source/assembly/locus,
                 selected model, sample and raw/true bounds. It subtracts the
                 supported model medians only for feature spans and targets;
                 evaluationseed is retained but no evaluation draws occur.
                 Valid centered rows retain the 69/22 numeric layout. The new
                 .centered.tsv ledgers retain EVERY original trial, its original
                 and centered bounds, applied medians and explicit invalid-span
                 or missing-feature holds. The ledger also records observed
                 model support, the support floor and the loaded family minimum
                 length. Crossed, below-minimum and out-of-context centers revert
                 to the original raw pair; unavailable context is not silently
                 labelled as known contig-edge clipping. A valid center with
                 missing NN features keeps its center without a further NN move.
                 Empty evaluation matrices are allowed
                 only when every frozen trial is explicitly held. Use the ledger
                 when comparing against the original raw baseline; the ordinary
                 scorer alone would measure its baseline in the centered frame.

Outputs: two full count/score tables and four compressed numeric matrices/radius,
conditioned-jitter accounting, complete original-locus ledger and population/stats.
Each worker selects a consensus once/locus and reuses it across all radii. Extra
features are candidate/selected-consensus length, whole-contig GC and candidate-span
Quantum identity. Each matrix has one compression thread in addition to t workers;
reserve enough CPUs for all concurrent writers. This command does not train networks,
call genomes, install resources or change annotations. Use -Xmx to set the heap.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.BuildNcrnaPositionalInputs "$@"
