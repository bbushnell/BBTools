#!/usr/bin/env bash
# Conditional matched-real-call endpoint comparison; no detection or training.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/evaluatencrnarealendpoints.sh config=literal_panel.config -Xmx2g -eoom

Required bindings:
 in=matched.witnesses.tsv gff=saved_off_calls.gff
 sources=source_roster.tsv sourceformat=native|legacy fastalist=ordered_fastas.txt
 familyroot=working_group ref=consensus.fa hbm=models.hbmc seeds=seeds.fa.gz
 family=cobalamin_rf00174|cobalamin_rf01482 (default: main member)
 mode=evaluate (default) or geometry
 inputs=panel_input_count witnesses=all_member_witness_count expected=member_count
 panel=panel_name out=fresh_prefix

Geometry mode writes out.geometry.tsv with exact assembly/call/model identity,
actual contig length, GFF fixed bits, gene-sense raw/truth ends and family length
bounds. It runs no alignment or inference. Do not pass offsets, tables, networks
or seedlabels in geometry mode. The cache is for conditional endpoint replay,
not caller reproduction or detection claims.
Legacy geometry additionally requires assemblyaudit=panel_sources.tsv: the
accepted outer audit maps exact FASTA paths to canonical assembly IDs, rather
than assuming legacy accession labels equal the IDs used for fold exclusion.

Gate-row mode: mode=gatevectors in=closed.geometry.tsv folds=assemblies.tsv
fold=0..4 (held internal assembly fold), or fold=-1 (complete outer panel).
Keep sources/sourceformat/fastalist/library/inputs/expected/panel bindings;
omit gff= and witnesses=. expected= is the entire cache count before selection.
Provide modeloffsets=fold_specific.offsets.tsv, newtables=six,paths and
newnets=start,stop,pairs with seedlabels=. Do not supply oldtables or oldnets.
Legacy sources still require the same assemblyaudit used for the cache.

This adapter selects real held calls before sequence extraction, rechecks FASTA
lengths and infers each available centered end once. Every selected call/end is
audited, including fixed, invalid-center, missing-feature and identical-raw holds.
Current-best/tied neural predictions may still produce an offset-only candidate.
The target is strictly closer to truth than RAW, not closer than the center.
Outputs: out.features.tsv, out.selection.tsv and, per sSEED.END, .vectors.tsv.gz
(35 inputs/one label), .rows.tsv.gz, .scores.tsv.gz and .summary.tsv. Compressed
outputs are BGZF. Raw FP32 scores are preserved for later cutoff/control replay.
Invalid single-end proposals stay labelled/audited; the subsequent evaluation
must apply actual single/joint geometry, reverting rejected ends to RAW.
No gate fit, cutoff selection, caller or new biological matching is performed.

Additional bindings for evaluate mode:
 modeloffsets=training_only.offsets.tsv
 oldtables=six,paths newtables=six,paths
 oldnets=start,stop,pairs newnets=start,stop,pairs seedlabels=58,59,60

Tables are ordered start probe0/1/2, stop probe0/1/2. Networks are ordered
start/stop within each declared seed. A centered network may be 'none'; that
seed/end is omitted from BOTH network comparisons, never filled by a checkpoint.
All available ends receive SINGLE_END_OTHER_AT_BASE rows, with the other end
held at raw (uncentered) or guarded center (centered). Only complete pairs also
receive PAIRED_JOINT rows. Single-end rows leave both-end accuracy empty. Both
scopes reuse the same saved inference outputs; never pool their denominators.
Both supported members use their ordinary 109-base minimum. Neural evaluation
retains radius ten, 69 inputs and 22 outputs.
Margins 0 and 0.1 reuse each available network evaluation. Outputs retain raw
FP32 feature/score bits, original real-call identities and every final coordinate.

Native sources use input_ordinal0/assembly_id/input_id/fasta/canonical_new_gff.
Legacy sources use the eleven-column NcrnaCalibrationPrep input roster; the
ordered FASTA list binds saved GFF input ordinals without filename inference.
Full winning model names must agree between witnesses, GFF and loaded library.
Fixed ends come from saved partial attributes, never merely from edge contact.
Invalid centers retain the raw pair and receive no centered inference. Missing
features retain the applicable base end. Native full-head confidence selection
and joint length/bounds holds apply. All original matched calls remain scored.

This is a CONDITIONAL endpoint-ability assay, not caller reproduction: saved
calls lack trim-success and NN-invocation provenance, and path competition and
detection are not rerun. The 1,000-genome rung is a less-rough comparison; the
legacy panel is regression evidence. Neither establishes adoption.

Finite controls: put 'test' first, followed by the same resource bindings and
out=fresh_directory. They construct tiny synthetic source-qualified fixtures.
HELP
 exit 0
fi
target=prok.EvaluateNcrnaRealEndpoints
if [[ $1 == test ]];then
 target=test.prok.EvaluateNcrnaRealEndpointsTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
