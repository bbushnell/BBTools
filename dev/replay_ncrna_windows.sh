#!/bin/bash
# Native saved-window acceptance development tool.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/replay_ncrna_windows.sh in=capture_prefix caller=exact_flat_baseline.config out=fresh_prefix expectedinputs=N expectedwindows=N expectedscheduled=N expectedmodels=N t=16 -Xmx2g -eoom
Replays EVERY saved RF00174 candidate, including genomic background, through the
real selected unlimited verifier at identity/borderline pairs60%,62%and64%.
Uses the actual caller configuration's resources and resolved controls.
Requires the original64% ordinary single-pass,109-base,two-distinct-key baseline.
Preserves original input order, zero-window inputs, physical contig bounds and
both coordinate systems. New .context.tsv/.context.fa.gz sidecars supply actual
refinement flanks. Legacy captures first use saved overlaps, then read missing
flanks from the exact original input, at most one streamed read per input.
The64% control must reproduce EACH original accepted proposal and alignment count.
Reports out.060/.062/.064 .results.tsv, .inputs.tsv and .replay.tsv.
Alignment totals count actual calls. Source-I/O totals are shared across all
three cuts and MUST NOT be added across cut summaries.
No final path is rerun and no final-call or precision result is invented.
Use summarize_ncrna_roster_baseline.sh acceptanceonly=t for original-gene
attribution, required64% global gene-set parity, coverage/truncation subsets
and unadjudicated background acceptances.
For G11's separate verifier-null control, use shufflewindows=10000 shuffleseed=1701.
This samples saved candidates without replacement and permutes their bases.
Discovery seeds are explicitly BYPASSED. Outside-window real refinement flanks
and original physical contig ends stay unchanged. Dinucleotide shuffle is absent.
The real patience-library all-model verifier runs at .60/.62/.64/.66/.68/.70,
with the109-base span floor on and off. Exact primary identities are reused
across settings; actual refinement/edge work and proposals are never cached.
Outputs .shuffle.windows.tsv and .shuffle.summary.tsv. Best PRIMARY identities
are reported separately from actual refined acceptance fractions over ALL sampled
windows. Logical requests and physical kernels are separate; no final path runs.
shufflecache=f provides the finite uncached correctness control. Ordinary replay
outputs/parity are not claimed or generated in this explicitly sampled null mode.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
	-cp "${DIR}current/" prok.ReplayNcrnaWindows "$@"
