#!/bin/bash
# Preserve every saved covariance-model hit and project its genomic relationships.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 6, 2026

Usage: dev/ncrnacmhitprojector.sh in=historical.tsv.gz contexts=crosswalk.tsv tables=searches.tsv out=new_directory

Projects every Infernal 1.1.5 cmsearch hit into genomic coordinates. Reads existing
tables only: no alignment, extraction, admission or SINGLE/TANDEM classification.
One metadata worker performs the joins; native I/O handles compression.

in=               Complete 21-column NcrnaAnnotatedCorpus locus table.
                  Contiguous intervals only; multipart annotations fail loudly.
contexts=         Ten-column sense-context map, one unique context ID per row:
                  context_id original_id genome seqid strand start1 stop1
                  contig_length lflank rflank (tab-separated header).
tables=           Search manifest with tab-separated columns:
                  table_id model_id model_sha80 clen ga_bits arm tblout
                  Model sha80 identifies the pinned model; the tool does not hash
                  or verify model bytes. Each table/path occurs once. All raw rows
                  are kept, including below-threshold, truncated and '?' hits.
out=              Fresh directory under an existing parent; overwrite is refused.
expectedloci=     Optional exact historical-annotation count, including unqueried rows.
expectedcontexts= Optional exact queried-context count, including zero-hit queries.
expectedhits=     Optional exact raw-hit count across all declared tables.

Coordinates are one-based inclusive. Flanks are actual sense-oriented lengths.
Genomic strand composes local hit strand with context orientation. Exact assembly,
contig, strand, model/version and genomic bounds define a physical hit interval;
distinct bounds remain distinct. No scores or E-values are aggregated.

Outputs:
events.tsv        Every raw event, source table/line, originating context and hit fields.
hits.tsv          Physical intervals with observation and annotation-overlap counts.
overlaps.tsv      All hit/annotation overlaps, including opposite-strand evidence.
queries.tsv       Every queried context and its event count; zero is not absence.
neighbors.tsv     Same-strand historical predecessor/successor in genomic start,
                  stop, ID order, plus all overlapping historical annotations.
                  Unqueried neighbors are explicitly labeled. Signed gaps are
                  successor.start-origin.stop-1 or origin.start-predecessor.stop-1;
                  negative means overlap. OVERLAP rows use NA for signed gap.
searches.tsv      Declared search/model identities and original table paths.
summary.tsv       Annotation, context, table, event and physical-interval counts.

Physical IDs are stable under input row reordering for a fixed set of intervals;
adding/removing intervals may renumber them. Their full tuple is in hits.tsv.
Standard BBTools config=, outstream= and JVM options are supported.
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
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaCmHitProjector "$@"
