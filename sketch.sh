#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 28, 2020

Description:  Creates one or more sketches from a fasta file,
optionally annotated with taxonomic information.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  sketch.sh in=<fasta file> out=<sketch file>

Standard parameters:
in=<file>           A fasta file containing one or more sequences.
out=<file>          Output filename.  If multiple files are desired it must
                    contain the # symbol.
blacklist=<file>    Ignore keys in this sketch file.  Additionally, there are
                    built-in blacklists that can be specified:
                       nt:      Blacklist for nt
                       refseq:  Blacklist for Refseq
                       silva:   Blacklist for Silva
                       img:     Blacklist for IMG
files=1             Number of output sketch files to produce, for parallel 
                    loading.  Independent of the number of sketches produced; 
                    sketches will be randomly distributed between files.
k=32,24             Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used, e.g. k=32,24
                    Query and reference k must match.
rcomp=t             Look at reverse-complement kmers also.
amino=f             Use amino acid mode.  Input should be amino acids.
translate=f         Call genes and translate to proteins.  Input should be
                    nucleotides.  Designed for prokaryotes.
mode=single         Possible modes:
                       single: Write one sketch.
                       sequence: Write one sketch per sequence.
                       taxa: Write one sketch per taxonomic unit.
                          Requires more memory, and taxonomic annotation.
                       img: Write one sketch per IMG id.
delta=t             Delta-compress sketches.
a48=t               Encode sketches as ASCII-48 rather than hex.
depth=f             Track the number of times kmers appear.
                    Required for the depth2 field in comparisons.
entropy=0.66        Ignore sequence with entropy below this value.
ssu=t               Scan for and retain full-length SSU sequence.

Size parameters:
size=10000          Desired size of sketches (if not using autosize).
maxfraction=0.01    (mgf) Max fraction of genomic kmers to use.
minsize=100         Do not generate sketches for genomes smaller than this.
autosize=t          Use flexible sizing instead of fixed-length.  This is
                    nonlinear; a human sketch is only ~6x a bacterial sketch.
sizemult=1          Multiply the autosized size of sketches by this factor.
                    Normally a bacterial-size genome will get a sketch size
                    of around 10000; if autosizefactor=2, it would be ~20000.
density=            If this flag is set (to a number between 0 and 1),
                    autosize and sizemult are ignored, and this fraction of
                    genomic kmers are used.  For example, at density=0.001,
                    a 4.5Mbp bacteria will get a 4500-kmer sketch.

Metadata parameters (optional; intended for single-sketch mode):
taxid=-1            Set the NCBI taxid.
imgid=-1            Set the IMG id.
spid=-1             Set the JGI sequencing project id.
name=               Set the name (taxname).
name0=              Set name0 (normally the first sequence header).
fname=              Set fname (normally the file name).
meta_=              Set an arbitrary metadata field.
                    For example, meta_Month=March.

Taxonomy-specific parameters:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.
imgdump=            Specify an IMG dump file containing NCBI taxIDs,
                    for IMG mode.
taxlevel=subspecies Taxa hits below this rank will be promoted and merged
                    with others.
prefilter=f         For huge datasets full of junk like nt, this flag
                    will save memory by ignoring taxa smaller than minsize.
                    Requires taxonomic information (tree and gi).
tossjunk=t          For taxa mode, discard taxonomically uninformative
                    sequences.  This includes sequences with no taxid,
                    with a tax level NO_RANK, of parent taxid of LIFE.
silva=f             Parse headers using Silva or semicolon-delimited syntax.

Ribosomal parameters, which allow SSU sequences to be attached to sketches:
processSSU=t        Run gene-calling to detect ribosomal SSU sequences.
16Sfile=<file>      Optional file of 16S sequences, annotated with TaxIDs.
18Sfile=<file>      Optional file of 18S sequences, annotated with TaxIDs.
preferSSUMap=f      Prefer file SSUs over called SSUs.
preferSSUMapEuks=t  Prefer file SSUs over called SSUs for Eukaryotes.
SSUMapOnly=f        Only use file SSUs.
SSUMapOnlyEuks=f    Only use file SSUs for Eukaryotes.  This prevents
                    associating an organism with its mitochondrial or
                    chloroplast 16S/18S, which is otherwise a problem.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx4g"
z2="-Xms4g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sketch() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP sketch.SketchMaker $@"
	echo $CMD >&2
	eval $CMD
}

sketch "$@"
