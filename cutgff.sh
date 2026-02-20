#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 15, 2019

Description:  Cuts out features defined by a gff file, and writes them
to a new fasta.  Features are output in their sense strand.

Usage:  cutgff.sh in=<fna file> gff=<gff file> out=<fna file>

in= is optional, and gff filenames will be automaitically assumed based on
the fasta name if not specified.  This allows running on multiple files:

cutgff.sh types=rRNA out=16S.fa minlen=1440 maxlen=1620 attributes=16S bacteria/*.fna.gz


File Parameters:
in=<file>           Input FNA (fasta) file.
gff=<file>          Input GFF file (optional).
out=<file>          Output FNA file.

Other Parameters:
types=CDS           Types of features to cut.
invert=false        Invert selection: rather outputting the features,
                    mask them with Ns in the original sequences.
attributes=         A comma-delimited list of strings.  If present, one of
                    these strings must be in the gff line attributes.
bannedattributes=   A comma-delimited list of banned strings.
banpartial=t        Ignore lines with 'partial=true' in attributes.
minlen=1            Ignore lines shorter than this.
maxlen=2147483647   Ignore lines longer than this.
renamebytaxid=f     Rename sequences with their taxID.  Input sequences
                    must be named appropriately, e.g. in NCBI format.
taxmode=accession   Valid modes are:
                       accession: Sequence names must start with an accession.
                       gi:        Seqence names must start with gi|number
                       taxid:     Sequence names must start with tid|number
                       header:    Best effort for various header formats.
requirepresent=t    Crash if a taxID cannot be found for a sequence.
oneperfile=f        Only output one sequence per file.
align=f             Align ribosomal sequences to consensus (if available);
                    discard those with low identity, and flip those
                    annotated on the wrong strand.
maxns=-1            If non-negative, ignore features with more than this many
                    undefined bases (Ns or IUPAC symbols).
maxnfraction=-1.0   If non-negative, ignore features with more than this
                    fraction of undefined bases (Ns or IUPAC symbols).
                    Should be 0.0 to 1.0.
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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP gff.CutGff $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
