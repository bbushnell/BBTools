#!/bin/bash

usage(){
echo "
Written by G11
Last modified September 2026

Description:  Runs NcrnaGffFamilyGrader -- coordinate-based, family-aware, one-to-one
6S (RF00013/RF01685) GFF grader. Grades a real callgenes.sh-produced GFF against a
coordinate-level truth manifest using maximum-cardinality bipartite matching (Kuhn's
algorithm), not the sealed NcrnaCombinedGradingDriver's per-record whole-count model --
this tool handles a real genome carrying multiple true loci on one contig. Writes
.rawcalls.tsv (a full per-locus audit trail: subtype outcome, pooled-union outcome,
cross-label diagnostic status, and stable call<->truth match IDs -- every aggregate in
.summary.tsv is reconstructable from this file alone) and .summary.tsv (per-family and
union TP/FN/FP/precision/recall, duplicate-query and negative-seqid counts, and the
cross-label diagnostic).

Usage:  ncrnagfffamilygrader.sh gff=<calls.gff> truth=<truth.tsv> [negseqids=<file>] out=<prefix>

Parameters:
gff=            A real callgenes.sh-produced GFF to grade.
truth=          Coordinate-level truth manifest: optional leading '#' comment lines,
                a header row, then EXACTLY 5 tab-delimited fields per row --
                seqid  family  start  stop  strand (1-based inclusive; family is
                exactly RF00013 or RF01685).
negseqids=      Optional: one seqid per line, known to carry zero true 6S loci --
                calls there are unmatched-by-construction straight FPs.
out=            Output prefix; writes <out>.rawcalls.tsv and <out>.summary.tsv.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

resolveSymlinks(){
	if [ -d "$DIR""current/" ]; then
		CP="$DIR""current/"
	elif [ -d "$DIR""../current/" ]; then
		CP="$DIR""../current/"
	fi
}

setEnv(){
	#DIR must be reassigned to the BBTools ROOT before sourcing javasetup.sh -- see the fix note in
	#testtrnakmerindex.sh (javasetup.sh reads DIR itself to find memdetect.sh).
	DIR="$DIR""../"
	. "$DIR""javasetup.sh"

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.NcrnaGffFamilyGrader $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
