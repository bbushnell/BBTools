#!/bin/bash
# Resolve primary ncRNA occurrences from existing covariance-model observations.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 7, 2026

Usage: ncrnaoccurrences.sh in=census primary=RF00504 out=new_directory

Resolves primary-model physical intervals into occurrences. An interval collapses
only when a parent contains its genomic bounds AND every child model-span event
has an actual containing parent event. Incomparable maximal parents are an error.
Other models supply overlap evidence only, never new or merged occurrences.
There is no sequence search, score/coverage threshold, or completeness filter.

in=          Accepted census directory: originals.tsv, searches.tsv, hits.tsv,
             events.tsv and overlaps.tsv from ncrnaobservationcensus.sh.
primary=     Required exact primary model ID, with one version in searches.tsv.
out=         Fresh directory below an existing parent; no overwrite or append.

Preserves original representative bounds, partial observations, every physical
interval/event link and every original annotation, including zero-hit originals.
Coordinates are one-based inclusive. Full genome/assembly/contig identity and
strand separate occurrences. Input tables and raw search evidence remain required
provenance; this tool writes joins to them, not replacement raw evidence.

One metadata worker. t= controls native I/O/compression settings, not independent
containment workers. Heap defaults to 1 GiB; all five metadata tables are retained.
Standard BBTools config=, outstream= and JVM options are supported.
See docs/guides/NcrnaOccurrencesGuide.txt for schemas and limitations.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaOccurrences "$@"
