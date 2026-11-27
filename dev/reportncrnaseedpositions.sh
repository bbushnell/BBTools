#!/usr/bin/env bash
# Empirical seed-edge caps and reference-only window comparisons.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 cat <<'HELP'
Usage: dev/reportncrnaseedpositions.sh in=member.fa seeds=forward17mers.fa.gz member=RF00174 expectedgenes=91072 expectedkeys=1975 orientation=sense out=fresh_prefix t=1 -Xmx1g

Serving-bound usage: add mode=bounds to the position-mode arguments above.
Optional boundpercentile=90|95|98|99 selects the nearest rank; default99.
The explicit value is recorded in every observed bound row and must match the
caller's keypercentile setting. This option is rejected outside bounds mode.
One counting pass over all admitted original genes writes .bounds.tsv and
.lengths.tsv only. Each key end uses the selected nearest-rank percentile
(default p99) plus10bases,capped by
the maximum length of the FULL member cohort minus17. Repeated occurrences
contribute separately; the fixed forward key axis and unobserved NA rows remain.
These margins are not the default position-mode p99.9 diagnostics. No second
window pass,caller,alignment,seed selection or model fitting runs in this mode.

Length-only usage: dev/reportncrnaseedpositions.sh mode=lengths in=genes.fa member=family out=fresh_prefix t=1 -Xmx1g
Reads complete FASTA records once and writes only .lengths.tsv: observed record
count, min/median/p99/p99.9/max and p99.9-minus17/max-minus17 padding caps.
expectedgenes=N is optional in this mode; when supplied it must match exactly.
Input may be stdin.fa, including a named FASTA streamed from an archive.
No seeds are needed; duplicate record names and sequence multiplicity are
preserved. These are reference-record lengths, not caller accuracy.

Streams already sense-oriented reference genes twice; does not reverse-complement
minus-strand annotations. No whole genomes, alignments, caller or seed selection.
Counts all occurrences and distinct references separately; retains zero keys.
Distances are from the key START to the gene5prime end and from the key END to
the gene3prime end. Caps use nearest-rankp99.9; n<1000 therefore uses the maximum.
The gene-length cap is nearest-rankp99.9 minus17, bounded below by zero.

A per-hit cap window is[keyStart-leftCap,keyEnd+rightCap], inclusive. Repeated
occurrences of one key form one envelope, then distinct-key envelopes combine
by either union or intersection within the known reference gene. No genomic
clipping, neighboring hits, caller clustering, claims, or length gates are run.
Zero-key references and one-key versus two-plus-key strata remain explicit.
Empty intersections have measured length0; absent windows and zero-key caps
are NA. Current native blanket geometry is[center-padding,center+padding+17],
so padding450 means918bases before clipping; this same formula tests the new
member length cap. Outputs are empirical in-corpus diagnostics, not caller recall.

Writes .keys.tsv, .lengths.tsv, .windows.tsv, .window_histogram.tsv and
.references.tsv. Quantiles use rank ceil(p*n), one-based, without interpolation.
Per-reference and histogram outputs may be large and should remain on the cluster.
One counting worker uses native streaming I/O; t=1 and complete unsampled inputs
are required. This diagnostic does not write a serving SeedOffsetTable resource.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ReportNcrnaSeedPositions "$@"
