#!/usr/bin/env bash
# Fixed-cohort coordinate replay from existing winner-qualified witness tables.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnamodelendpointoffsets.sh in=train.witnesses.tsv \
 refs=main.fa,second.fa,third.fa members=RF00174,RF01482,RF01689 family=cobalamin \
 evaluation=new1000.witnesses.tsv,old300.witnesses.tsv panels=new1000,old300 \
 bindings=training.sources.tsv exclude=outer.panel_sources.tsv \
 expectedmodels=412 expected=19995,17863,2036,380 thresholds=5,10,20 out=report

Learns each winning model's median signed error using only training assemblies
outside both outer panels. in= may name comma-separated disjoint witness files.
Expected counts: training before exclusions, retained training, then each panel.
All catalog models, including zero-support models, enter models/support tables.
Quantiles use nearest rank, with the lower integer median when n is even.

Replays held-out errors after subtracting training medians. Unsupported models
hold. Invalid one-based/joint intervals hold both ends and are counted.
Output: .models.tsv, .support.tsv, .population.tsv, .single.tsv, .both.tsv.
Tables include exact/within3/within10 counts and percentages, RMS before/after/
signed change, supported/shifted/geometry-held counts and better/worse/tied.
Both-end tolerances require both ends of one locus. Joint RMS is pooled over
two ends per locus; joint better/worse compares their summed squared error.
No rematching, detection claim, contig-length clipping, gene-path competition,
network training, calling, or resource adoption occurs. Outputs must be fresh.

Optional breakdown=RF00174 breakdowncounts=1718,286 thresholds=10 writes ONLY
.breakdown.tsv, partitioning that member's fixed matches by applied median size
and training IQR: 0,1-2,3-5,6-10,above10 (missing IQR reported separately).
Single-end bins use that end; both-end bins use max(abs(median5),abs(median3))
and max(IQR5,IQR3). Unsupported and invalid spans apply zero shift; support and
geometry remain explicit. This is descriptive and never selects a shifting rule.

Optional centeredhistograms=t thresholds=10 writes ONLY one signed residual
histogram per member, .offsets.tsv and .training.tsv. It learns medians before
a second small training-witness pass, excluding the same outer assemblies and
preserving joint geometry holds. The compact resource keeps every ordered model,
its support and applied signed medians; unsupported models have zero shifts.
These histograms feed jitter sampling and weighted-background regeneration;
they contain no evaluation matches. They do not fit or serve a network.
For a minimum-length-aware residual export, add minlengths=L1,L2,L3 in the
same order as members=. Bind each value to that member's loaded family minimum.
Below-minimum centered spans retain BOTH original errors, with coordinate holds,
minimum-length holds and already-short raw spans counted separately in training.tsv.
Omitting minlengths preserves the historical coordinate-only export byte-for-byte.
Whole-contig lengths are unavailable in witnesses; no upper-contig-bound claim
is made. This option applies only to centeredhistograms, not historical replay.

All-data product offsets and residuals:
refit=t centeredhistograms=t thresholds=10 expected=TOTAL,TOTAL
Omit exclude=, evaluation= and panels=. Every source-bound saved match teaches,
including former held-out assemblies. Also writes models/support/population
tables for stability review, explicitly labeled no-holdout/no-evaluation;
no single/both held-out replay tables are produced. Source bindings stay required.

Saved-table stability only:
dev/ncrnamodelendpointoffsets.sh compare=old.models.tsv,new.models.tsv \
 members=RF00174,RF01482,RF01689 family=cobalamin expectedmodels=412 \
 thresholds=10 out=comparison
Reads only the two complete model tables on the same ordered axis. New support
must include old support. Writes .stability.tsv with signed new-minus-old median
changes and .stability_summary.tsv with absolute-change bins among models meeting
minimum support in BOTH populations. Newly supported models are counted separately;
missing old medians and their changes remain NA. No witnesses or panels are read.

Finite controls: dev/ncrnamodelendpointoffsets.sh test out=fresh_directory
HELP
 exit 0
fi
target=prok.NcrnaModelEndpointOffsets
if [[ $1 == test ]];then
 target=test.prok.NcrnaModelEndpointOffsetsTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
