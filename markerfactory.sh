#!/bin/bash

usage(){
echo "
Written by Eru
Last modified July 31, 2026

Description:  Builds per-domain single-copy marker sets from a manifest of
per-genome protein FASTAs.  Pools all proteins, clusters them into families
(greedy identity-threshold clustering, CD-HIT / linclust shape), then per
domain computes each family's prevalence and full copy-number distribution
(how many genomes carry it exactly 0/1/2/3/4+ times).  A family is selected as
a single-copy marker when it is carried EXACTLY ONCE in at least 'threshold'
fraction of that domain's genomes.  Writes the marker set as a TSV plus a
provenance sidecar (out.prov).

Usage:  markerfactory.sh manifest=<genomes.tsv> out=<markers.tsv>

Manifest format:
Tab-separated text, one genome per line, columns:
genome_id <tab> domain <tab> lineage <tab> fasta_path
Blank lines and lines starting with # are ignored.

Required parameters:
manifest=<file> Genome manifest (see format above).

Optional parameters (and their defaults):
out=stdout      Output TSV.  A .prov provenance sidecar is written beside it.
threshold=0.97  Min fraction of a domain's genomes carrying a family EXACTLY
                ONCE to select it as a single-copy marker.
minid=90        Clustering min percent identity (0.9 also accepted).
mincov=0.8      Clustering min aligned fraction of both member and rep.
k=5             Clustering seed k-mer length.
reduced=f       Use amino8 reduced-alphabet seeds (more sensitive).
version=v1      Marker-set version id.
timestamp=NA    Build timestamp recorded in provenance (not read from clock).
taxonomy=NA     Taxonomy snapshot id recorded in provenance.
ow=t            (overwrite) Overwrite existing output.

Output columns (tab-separated):
domain marker_set_version family_id representative prevalence n_genomes
copies_0 copies_1 copies_2 copies_3 copies_4plus fraction_exactly_once
selected_single_copy

Notes:
One row per family per domain; families with zero prevalence in a domain are
omitted.  Family ids are stable within a run only (cross-run stable ids and
the lineage hierarchy with nearest-ancestor fallback are not implemented).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"

}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.MarkerFactoryCLI $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
