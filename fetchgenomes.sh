#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 10, 2026

Description:  Selects a phylogenetically diverse set of RefSeq genomes from NCBI
assembly_summary.txt files and writes a hardened download script (curl with retries,
gzip integrity checks, fna+gff pairing, _tid_ header renaming). Replacement for
fetchproks.sh; no FTP crawling — selection happens locally from the summary table.

Get summaries first, e.g.:
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
curl -o archaea_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt

Usage:  fetchgenomes.sh summary=<file[,file2]> out=<script.sh> [options]

Parameters:
summary=<file>      Comma-list of assembly_summary.txt files (required).
out=<file>          Output download script (required).
tree=f              Load TaxTree for lineage-aware quotas (default f); tree=t, tree=auto,
                   or tree=<path> enables it.  treefile=<path> selects the source only.
usetree=f           Alias for tree=f.
maxperspecies=1     Assemblies kept per species (best-ranked first).
maxpergenus=2       Species kept per genus (0=unlimited).
maxperfamily=0      Species kept per family (0=unlimited).
skip=<file>         One taxid per line to exclude (e.g. an existing collection).
allowexcluded=f     Keep assemblies flagged excluded_from_refseq.
rename=t            Add _tid_<taxid> suffix to sequence headers in the script.
minsize=0           Minimum genome_size when the column is present.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.FetchGenomes $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
