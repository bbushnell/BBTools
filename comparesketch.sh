#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 7, 2020

Description:  Compares query sketches to others, and prints their kmer identity.
The input can be sketches made by sketch.sh, or fasta/fastq files.
It's recommended to first sketch references with sketch.sh for large files,
or when taxonomic information is desired.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  comparesketch.sh in=<file,file,file...> ref=<file,file,file...>
Alternative:  comparesketch.sh in=<file,file,file...> file file file
Alternative:  comparesketch.sh in=<file,file,file...> *.sketch
Alternative:  comparesketch.sh alltoall *.sketch.gz

File parameters:
in=<file,file...>   Sketches or fasta files to compare.
out=stdout          Comparison output.  Can be set to a file instead.
outsketch=<file>    Optionally write sketch files generated from the input.
ref=<file,file...>  List of sketches to compare against.  Files given without
                    a prefix (ref=) will be treated as references,
                    so you can use *.sketch without ref=.
                    You can also do ref=nt#.sketch to load all numbered files
                    fitting that pattern.
                    On NERSC, you can use these abbreviations (e.g. ref=nt):
                       nt:      nt sketches
                       refseq:  Refseq sketches
                       silva:   Silva sketches
                       img:     IMG sketches
                       mito:    RefSeq mitochondrial sketches
                       fungi:   RefSeq fungi sketches
                       protein: RefSeq prokaryotic amino acid sketches
                    Using an abbreviation automatically sets the blacklist, 
                    and k.  If the reference is in amino space, the query
                    also be in amino acid space with the flag amino added.
                    If the query is in nucleotide space, use the flag
                    'translate', but this will only work for prokaryotes.

Blacklist and Whitelist parameters:
blacklist=<file>    Ignore keys in this sketch file.  Additionally, there are
                    built-in blacklists that can be specified:
                       nt:      Blacklist for nt
                       refseq:  Blacklist for Refseq
                       silva:   Blacklist for Silva
                       img:     Blacklist for IMG
whitelist=f         Ignore keys that are not in the index.  Requires index=t.

Sketch-making parameters:
mode=perfile        Possible modes, for sequence input:
                       single: Generate one sketch.
                       sequence: Generate one sketch per sequence.
                       perfile: Generate one sketch per file.
sketchonly=f        Don't run comparisons, just write the output sketch file.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
                    Dual kmers are fastest if the shorter is a multiple 
                    of 4.  Query and reference k must match.
samplerate=1        Set to a lower value to sample a fraction of input reads.
                    For raw reads (rather than an assembly), 1-3x coverage
                    gives best results, by reducing error kmers.  Somewhat
                    higher is better for high-error-rate data like PacBio.
minkeycount=1       Ignore kmers that occur fewer times than this.  Values
                    over 1 can be used with raw reads to avoid error kmers.
minprob=0.0001      Ignore kmers below this probability of correctness.
minqual=0           Ignore kmers spanning bases below this quality.
entropy=0.66        Ignore sequence with entropy below this value.
merge=f             Merge paired reads prior to sketching.
amino=f             Use amino acid mode.  Input should be amino acids.
translate=f         Call genes and translate to proteins.  Input should be
                    nucleotides.  Designed for prokaryotes.
sixframes=f         Translate all 6 frames instead of predicting genes.
ssu=t               Scan for and retain full-length SSU sequence.
printssusequence=f  Print the query SSU sequence (JSON mode only).

Size parameters:
size=10000          Desired size of sketches (if not using autosize).
mgf=0.01            (maxfraction) Max fraction of genomic kmers to use.
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
sketchheapfactor=4  If minkeycount>1, temporarily track this many kmers until
                    counts are known and low-count kmers are discarded.

Sketch comparing parameters:
threads=auto        Use this many threads for comparison.
index=auto          Index the sketches for much faster searching.
                    Requires more memory and adds startup time.
                    Recommended true for many query sketches, false for few.
prealloc=f          Preallocate the index for greater efficiency.
                    Can be set to a number between 0 and 1 to determine how 
                    much of total memory should be used.
alltoall            (ata) Compare all refs to all.  Must be sketches.
compareself=f       In all-to-all mode, compare a sketch to itself.

Taxonomy-related parameters:
tree=<file>         Specify a TaxTree file.  On Genepool, use tree=auto.
                    Only necessary for use with printtaxa and level.
                    Assumes comparisons are done against reference sketches
                    with known taxonomy information.
level=2             Only report the best record per taxa at this level.
                    Either level names or numbers may be used.
                        0: disabled
                        1: subspecies
                        2: species
                        3: genus
                       ...etc
include=            Restrict output to organisms in these clades.
                    May be a comma-delimited list of names or NCBI TaxIDs.
includelevel=0      Promote the include list to this taxonomic level.
                    For example, include=h.sapiens includelevel=phylum
                    would only include organisms in the same phylum as human.
includestring=      Only report records whose name contains this string.
exclude=            Ignore organisms in these clades.
                    May be a comma-delimited list of names or NCBI TaxIDs.
excludelevel=0      Promote the exclude list to this taxonomic level.
                    For example, exclude=h.sapiens excludelevel=phylum
                    would exclude all organisms in the same phylum as human.
excludestring=      Do not records whose name contains this string.
minlevel=           Use this to restrict comparisons to distantly-related
                    organisms.  Intended for finding misclassified organisms
                    using all-to-all comparisons.  minlevel=order would only
                    report hits between organisms related at the order level
                    or higher, not between same species or genus.
banunclassified=f   Ignore organisms descending from nodes like 
                    'unclassified Bacteria'
banvirus=f          Ignore viruses.
requiressu=f        Ignore records without SSUs.
minrefsize=0        Ignore ref sketches smaller than this (unique kmers).
minrefsizebases=0   Ignore ref sketches smaller than this (total base pairs).

Output format parameters:
format=2            2: Default format with, per query, one query header line;
                       one column header line; and one reference line per hit.
                    3: One line per hit, with columns query, reference, ANI,
                       and sizeRatio. Useful for all-to-all comparisons.
                    4: JSON (format=json also works).
                    5: Constellation (format=constellation also works).
usetaxidname=f      For format 3, print the taxID in the name column.
usetaxname          for format 3, print the taxonomic name in the name column.
useimgname          For format 3, print the img ID in the name column.

Output column parameters (for format=2):
printall=f          Enable all output columns.
printani=t          (ani) Print average nucleotide identity estimate.
completeness=t      Genome completeness estimate.
score=f             Score (used for sorting the output).
printmatches=t      Number of kmer matches to reference.
printlength=f       Number of kmers compared.
printtaxid=t        NCBI taxID.
printimg=f          IMG identifier (only for IMG data).
printgbases=f       Number of genomic bases.
printgkmers=f       Number of genomic kmers.
printgsize=t        Estimated number of unique genomic kmers.
printgseqs=t        Number of sequences (scaffolds/reads).
printtaxname=t      Name associated with this taxID.
printname0=f        (pn0) Original seqeuence name.
printfname=t        Query filename.
printtaxa=f         Full taxonomy of each record.
printcontam=t       Print contamination estimate, and factor contaminant kmers
                    into calculations.  Kmers are considered contaminant if
                    present in some ref sketch but not the current one.
printunique=t       Number of matches unique to this reference.
printunique2=f      Number of matches unique to this reference's taxa.
printunique3=f      Number of query kmers unique to this reference's taxa,
                    regardless of whether they are in this reference sketch.
printnohit=f        Number of kmers that don't hit anything.
printrefhits=f      Average number of ref sketches hit by shared kmers.
printgc=f           GC content.
printucontam=f      Contam hits that hit exactly one reference sketch.
printcontam2=f      Print contamination estimate using only kmer hits
                    to unrelated taxa.
contamlevel=species Taxonomic level to use for contam2/unique2/unique3.
NOTE: unique2/unique3/contam2/refhits require an index.

printdepth=f        (depth) Print average depth of sketch kmers; intended
                    for shotgun read input.
printdepth2=f       (depth2) Print depth compensating for genomic repeats.
                    Requires reference sketches to be generated with depth.
actualdepth=t       If this is false, the raw average count is printed.
                    If true, the raw average (observed depth) is converted 
                    to estimated actual depth (including uncovered areas).
printvolume=f       (volume) Product of average depth and matches.
printca=f           Print common ancestor, if query taxID is known.
printcal=f          Print common ancestor tax level, if query taxID is known.
recordsperlevel=0   If query TaxID is known, and this is positive, print this
                    many records per common ancestor level.

Sorting parameters:
sortbyscore=t       Default sort order is by score, a composite metric.
sortbydepth=f       Include depth as a factor in sort order.
sortbydepth2=f      Include depth2 as a factor in sort order.
sortbyvolume=f      Include volume as a factor in sort order.
sortbykid=f         Sort strictly by KID.
sortbyani=f         Sort strictly by ANI/AAI/WKID.
sortbyhits=f        Sort strictly by the number of kmer hits.

Other output parameters:
minhits=3           (hits) Only report records with at least this many hits.
minani=0            (ani) Only report records with at least this ANI (0-1).
minwkid=0.0001      (wkid) Only report records with at least this WKID (0-1).
anifromwkid=t       Calculate ani from wkid.  If false, use kid.
minbases=0          Ignore ref sketches of sequences shortert than this.
minsizeratio=0      Don't compare sketches if the smaller genome is less than
                    this fraction of the size of the larger.
records=20          Report at most this many best-matching records.
color=family        Color records at the family level.  color=f will disable.
                    Colors work in most terminals but may cause odd characters
                    to appear in text editors.  So, color defaults to f if 
                    writing to a file.  Requires the taxtree to be loaded.
intersect=f         Print sketch intersections.  delta=f is suggested.

Metadata parameters (optional, for the query sketch header):
taxid=-1            Set the NCBI taxid.
imgid=-1            Set the IMG id.
spid=-1             Set the JGI sequencing project id.
name=               Set the name (taxname).
name0=              Set name0 (normally the first sequence header).
fname=              Set fname (normally the file name).
meta_=              Set an arbitrary metadata field.
                    For example, meta_Month=March.

Other parameters:
requiredmeta=       (rmeta) Required optional metadata values.  For example:
                    rmeta=subunit:ssu,source:silva
bannedmeta=         (bmeta) Forbidden optional metadata values.


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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

comparesketch() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP sketch.CompareSketch $@"
#	echo $CMD >&2
	eval $CMD
}

comparesketch "$@"
