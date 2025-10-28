#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 10, 2024

Description:  Compares query sketches to reference sketches hosted on a 
remote server via the Internet.  The input can be sketches made by sketch.sh,
or fasta/fastq files from which SendSketch will generate sketches.  
Only sketches will sent, not sequences.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  
sendsketch.sh in=file

To change nucleotide servers, add the server name, e.g.:
sendsketch.sh in=file nt

For the protein server with nucleotide input:
sendsketch.sh in=file protein

For the protein server with amino input:
sendsketch.sh in=file amino protein


Standard parameters:
in=<file>       Sketch or fasta file to compare.
out=stdout      Comparison output.  Can be set to a file instead.
outsketch=      Optional, to write the sketch to a file.
local=f         For local files, have the server load the sketches.
                Allows use of whitelists; recommended for Silva.
                Local can only be used when the client and server access 
                the same filesystem - e.g., Genepool and Cori.
address=        Address of remote server.  Default address:
                https://refseq-sketch.jgi.doe.gov/sketch
                You can also specify these abbreviations:
                   nt:      nt server
                   refseq:  Refseq server
                   silva:   Silva server
                   protein: RefSeq prokaryotic amino acid sketches
                   img:     IMG server (Not Yet Available)
                   mito:    RefSeq mitochondrial server (NYA)
                   fungi:   RefSeq fungi sketches (NYA)
                Using an abbreviation automatically sets the address, 
                the blacklist, and k.
aws=f           Set aws=t to use the aws servers instead of NERSC.
                When, for example, NERSC (or the whole SF Bay area) is down.

Sketch-making parameters:
mode=single     Possible modes, for fasta input:
                   single: Generate one sketch per file.
                   sequence: Generate one sketch per sequence.
k=31            Kmer length, 1-32.  This is automatic and does not need to
                be set for JGI servers, only for locally-hosted servers.
samplerate=1    Set to a lower value to sample a fraction of input reads.
                For raw reads (rather than an assembly), 1-3x coverage
                gives best results, by reducing error kmers.  Somewhat
                higher is better for high-error-rate data like PacBio.
minkeycount=1   Ignore kmers that occur fewer times than this.  Values
                over 1 can be used with raw reads to avoid error kmers.
minprob=0.0001  Ignore kmers below this probability of correctness.
minqual=0       Ignore kmers spanning bases below this quality.
entropy=0.66    Ignore sequence with entropy below this value.
merge=f         Merge paired reads prior to sketching.
amino=f         Use amino acid mode.  Input should be amino acids.
translate=f     Call genes and translate to proteins.  Input should be
                nucleotides.  Designed for prokaryotes.
sixframes=f     Translate all 6 frames instead of predicting genes.
ssu=t           Scan for and retain full-length SSU sequence.
printssusequence=f  Print the query SSU sequence (JSON mode only).
refid=          Instead of a query file, specify a reference sketch by name
                or taxid; e.g. refid=h.sapiens or refid=9606.

Size parameters:
size=10000      Desired size of sketches (if not using autosize).
mgf=0.01        (maxfraction) Max fraction of genomic kmers to use.
minsize=100     Do not generate sketches for genomes smaller than this.
autosize=t      Use flexible sizing instead of fixed-length.  This is
                nonlinear; a human sketch is only ~6x a bacterial sketch.
sizemult=1      Multiply the autosized size of sketches by this factor.
                Normally a bacterial-size genome will get a sketch size
                of around 10000; if autosizefactor=2, it would be ~20000.
density=        If this flag is set (to a number between 0 and 1),
                autosize and sizemult are ignored, and this fraction of
                genomic kmers are used.  For example, at density=0.001,
                a 4.5Mbp bacteria will get a 4500-kmer sketch.
sketchheapfactor=4  If minkeycount>1, temporarily track this many kmers until
                counts are known and low-count kmers are discarded.

Taxonomy and filtering parameters:
level=2         Only report the best record per taxa at this level.
                Either level names or numbers may be used.
                    0: disabled
                    1: subspecies
                    2: species
                    3: genus
                   ...etc
include=        Restrict output to organisms in these clades.
                May be a comma-delimited list of names or NCBI TaxIDs.
includelevel=0  Promote the include list to this taxonomic level.
                For example, include=h.sapiens includelevel=phylum
                would only include organisms in the same phylum as human.
includestring=  Only report records whose name contains this string.
exclude=        Ignore organisms in these clades.
                May be a comma-delimited list of names or NCBI TaxIDs.
excludelevel=0  Promote the exclude list to this taxonomic level.
                For example, exclude=h.sapiens excludelevel=phylum
                would exclude all organisms in the same phylum as human.
excludestring=  Do not records whose name contains this string.
banunclassified=f   Ignore organisms descending from nodes like 
                    'unclassified Bacteria'
banvirus=f      Ignore viruses.
requiressu=f    Ignore records without SSUs.
minrefsize=0    Ignore ref sketches smaller than this (unique kmers).
minrefsizebases=0   Ignore ref sketches smaller than this (total base pairs).

Output format parameters:
format=2        2: Default format with, per query, one query header line;
                   one column header line; and one reference line per hit.
                3: One line per hit, with columns query, reference, ANI,
                   and sizeRatio.
                4: JSON (format=json also works).
                5: Constellation (format=constellation also works).
usetaxidname=f  For format 3, print the taxID in the name column.
usetaxname      for format 3, print the taxonomic name in the name column.
useimgname      For format 3, print the img ID in the name column.
d3=f            Output in JSON format, with a tree for visualization.

Output column parameters (for format=2):
printall=f      Enable all output columns.
printani=t      (ani) Print average nucleotide identity estimate.
completeness=t  Genome completeness estimate.
score=f         Score (used for sorting the output).
printmatches=t  Number of kmer matches to reference.
printlength=f   Number of kmers compared.
printtaxid=t    NCBI taxID.
printimg=f      IMG identifier (only for IMG data).
printgbases=f   Number of genomic bases.
printgkmers=f   Number of genomic kmers.
printgsize=t    Estimated number of unique genomic kmers.
printgseqs=t    Number of sequences (scaffolds/reads).
printtaxname=t  Name associated with this taxID.
printname0=f    (pn0) Original seqeuence name.
printqfname=t   Query filename.
printrfname=f   Reference filename.
printtaxa=f     Full taxonomy of each record.
printcontam=t   Print contamination estimate, and factor contaminant kmers
                into calculations.  Kmers are considered contaminant if
                present in some ref sketch but not the current one.
printunique=t   Number of matches unique to this reference.
printunique2=f  Number of matches unique to this reference's taxa.
printunique3=f  Number of query kmers unique to this reference's taxa,
                regardless of whether they are in this reference sketch.
printnohit=f    Number of kmers that don't hit anything.
printrefhits=f  Average number of ref sketches hit by shared kmers.
printgc=f       GC content.
printucontam=f  Contam hits that hit exactly one reference sketch.
printcontam2=f  Print contamination estimate using only kmer hits
                to unrelated taxa.
contamlevel=species Taxonomic level to use for contam2/unique2/unique3.
NOTE: unique2/unique3/contam2/refhits require an index.

printdepth=f    (depth) Print average depth of sketch kmers; intended
                for shotgun read input.
printdepth2=f   (depth2) Print depth compensating for genomic repeats.
                Requires reference sketches to be generated with depth.
actualdepth=t   If this is false, the raw average count is printed.
                If true, the raw average (observed depth) is converted 
                to estimated actual depth (including uncovered areas).
printvolume=f   (volume) Product of average depth and matches.
printca=f       Print common ancestor, if query taxID is known.
printcal=f      Print common ancestor tax level, if query taxID is known.
recordsperlevel=0   If query TaxID is known, and this is positive, print at
                    most this many records per common ancestor level.

Sorting parameters:
sortbyscore=t   Default sort order is by score.
sortbydepth=f   Include depth as a factor in sort order.
sortbydepth2=f  Include depth2 as a factor in sort order.
sortbyvolume=f  Include volume as a factor in sort order.
sortbykid=f     Sort strictly by KID.
sortbyani=f     Sort strictly by ANI/AAI/WKID.
sortbyhits=f    Sort strictly by the number of kmer hits.

Other output parameters:
minhits=3       (hits) Only report records with at least this many hits.
minani=0        (ani) Only report records with at least this ANI (0-1).
minwkid=0.0001  (wkid) Only report records with at least this WKID (0-1).
anifromwkid=t   Calculate ani from wkid.  If false, use kid.
minbases=0      Ignore ref sketches of sequences shortert than this.
minsizeratio=0  Don't compare sketches if the smaller genome is less than
                this fraction of the size of the larger.
records=20      Report at most this many best-matching records.
color=family    Color records at the family level.  color=f will disable.
                Colors work in most terminals but may cause odd characters
                to appear in text editors.  So, color defaults to f if 
                writing to a file.
intersect=f     Print sketch intersections.  delta=f is suggested.

Metadata parameters (optional, for the query sketch header):
taxid=-1        Set the NCBI taxid.
imgid=-1        Set the IMG id.
spid=-1         Set the sequencing project id (JGI-specific).
name=           Set the name (taxname).
name0=          Set name0 (normally the first sequence header).
fname=          Set fname (normally the file name).
meta_=          Set an arbitrary metadata field.
                For example, meta_Month=March.

Other parameters:
requiredmeta=   (rmeta) Required optional metadata values.  For example:
                rmeta=subunit:ssu,source:silva
bannedmeta=     (bmeta) Forbidden optional metadata values.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sendsketch() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP sketch.SendSketch $@"
#	echo $CMD >&2
	eval $CMD
}

sendsketch "$@"
