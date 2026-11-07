#!/bin/bash
# Describe saved ncRNA observations and actual searched-context coverage.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 6, 2026

Usage: ncrnaobservationcensus.sh in=projection originals=loci.tsv contexts=crosswalk.tsv queryfiles=queries.tsv out=new_directory

Describes all saved observations from NcrnaCmHitProjector: physical overlaps,
annotation relationships, observed neighbors and actual searched coverage.
No extraction, alignment, model selection, admission or biological labels.
One metadata worker; native byte I/O and primitive indices. All input populations
are complete. Output must be a fresh directory beneath an existing parent.

in=          Directory containing events.tsv, hits.tsv, overlaps.tsv, queries.tsv,
             neighbors.tsv and searches.tsv from ncrnacmhitprojector.sh.
originals=   Original 21-column annotated-corpus metadata, with all flags intact.
contexts=    Original ten-column sense-context crosswalk used by the projector.
queryfiles=  Tab-separated header table_id query_fasta; exactly one binding for
             each search. Existing FASTA paths must match the corresponding
             Infernal table's '# Target file:' record. Read once per distinct
             path. Query IDs must be ASCII, unique within a file and present
             in the context map; sequences must have the declared context length.
out=         Fresh directory. Existing outputs are never overwritten or appended.

Original evidence rows are copied unchanged. New tables retain every overlapping
physical pair, each original annotation including zero-hit ones, all annotation
child pairs, and nearest nonoverlapping hit neighbors on each side and strand.
All ties are retained. Grouping uses genome, assembly AND contig.

Coordinates are one-based inclusive. Pair orientation is lexical strand/model/
version then numeric start/stop within a common source. Gap=right.start-left.stop-1;
adjacent intervals have gap zero. Pair fractions use EACH interval's own length.

Shared-context evidence requires one actual context searched by the relevant
model(s). Unions of multiple contexts are separate evidence and retain every
uncovered segment. Zero observations and missing neighbors never become negative
biological labels. Raw model coverage, scores, E-values and truncation remain
per-event and are never averaged or used to select a winner.

See docs/guides/NcrnaObservationCensusGuide.txt for output schemas and joins.
Standard BBTools config=, outstream= and JVM options are supported.
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter problems.
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
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaObservationCensus "$@"
