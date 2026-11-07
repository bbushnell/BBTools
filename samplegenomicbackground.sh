#!/bin/bash
# Samples annotation-masked genomic windows; masking alone is not a negative label.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Raiden
Last modified November 5, 2026

Usage: samplegenomicbackground.sh mode=panel genomes=train.genomes.tsv
       groups=source.groups.tsv expectedgenomes=28809 panelcount=1000 seed=47
       out=panel

Usage: samplegenomicbackground.sh genomes=all.genomes.tsv panelids=panel.ids.txt
       in=train.fa.gz loci=corpus.loci.tsv.gz families=RF00059 arm=panel
       expectedgenomes=39979 expectedpanel=1000 genomecount=1000
       expectedlengths=58125 expectedloci=70671 count=10000 flank=50 seed=101
       stratify=phylumgc gcbins=10 nonoverlap=t out=background

Selects source panels or samples real genomic windows with empirical lengths.
Examples illustrate schemas/options; replace counts with your actual populations.
Outputs are PROVISIONAL backgrounds, not confirmed negatives. Independently screen
the entire population before interpreting acceptance rates as false positives.

Shared options:
mode=sample       sample or panel.
genomes=          Original eleven-field NcrnaCorpusSample genome roster.
expectedgenomes=  Required exact number of rows in genomes=.
out=             Fresh output prefix under an existing parent (alias outprefix=).
seed=            Required signed integer; fixed inputs/order policy give fixed draws.
t=1              One owner; other thread counts rejected.
ow=f             Existing outputs are never overwritten.
-Xmx2g           Heap; sampling retains one source genome at a time.

Panel selection (no genome sequence reads):
groups=          Six-field native census source-group ledger.
panelcount=1000  Maximum selected source groups. Input roster must contain whole
                 eligible groups, already restricted to the permitted split.
                 Exclude taxonomy-ambiguous groups; retain at least one per phylum.
                 Allocate remaining slots by remaining group capacity and choose
                 groups with a seeded shuffle. Representative source is lexical.
Outputs: .ids.txt, .genomes.tsv, .groups.tsv, .summary.tsv.

Window sampling:
panelids=        Exact unique source IDs, one per line; no header or provenance hash.
panel=           Legacy two-column sha80/absolute-FASTA ledger instead of panelids=.
arm=             panel (all supplied panel sources) or nonpanel (legacy alias exclusion).
genomecount=     Required number of selected source genomes; panel uses expectedpanel.
expectedpanel=   Required exact number of supplied panel entries.
in=              FASTA of empirical gene lengths (alias lengths=); IDs must be unique.
expectedlengths= Required number of empirical FASTA records.
loci=            Comma-separated original 21-field source-locus ledgers.
families=        Matching Rfam IDs, one per locus ledger, in the same order.
                 Default: RF00174,RF01482,RF01689.
expectedloci=    Required total rows across all supplied locus ledgers.
count=           Required number of windows; must be at least genomecount.
flank=           Required nonnegative padding for contexts; clips at contig edges.
stratify=none    phylumgc balances windows across occupied domain/phylum/GC strata,
                 then sources within each stratum; requires arm=panel.
gcbins=10        Positive bin count, only with phylumgc. Equal-width half-open GC
                 fractions; 100% goes into the last bin. Counts exclude ambiguous bases.
nonoverlap=f     t excludes each selected core from later draws in its source contig.
                 Contexts may overlap. Exhaustion fails without redrawing lengths.

Default outputs: .cores.fa.gz, .contexts.fa.gz, .shreds.tsv, .genomes.tsv,
.lengths.tsv, .summary.tsv. Stratification/nonoverlap additionally write the entire
.schedule.tsv before sampling and .window_metadata.tsv; stratification also writes
.strata.tsv with selected versus contributing source counts. Coordinates are
one-based inclusive. Failed runs retain partial, closed outputs and the fixed schedule.
Defaults retain historical cobalamin sampling behavior and output schemas.

Full schemas and methodology: docs/guides/SampleGenomicBackgroundGuide.txt
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
parseJavaArgs --xmx=2g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.CobalaminShredSampler "$@"
