#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 12, 2024

Description:  Performs high-speed alignment-free sequence quantification,
by counting the number of long kmers that match between a read and
a set of reference sequences.  Designed for RNA-seq with alternative splicing.
Please read bbmap/docs/guides/SealGuide.txt for more information.

Usage:  seal.sh in=<file> *.fa pattern=out_%.fq outu=unmapped.fq stats=stats.txt

Sequence quantification examples:
seal.sh in=<file> ref=<file> rpkm=rpkm.txt stats=stats.txt
or
seal.sh in=<file> ref=<file,file,file...> refstats=refstats.txt

Splitting examples:
seal.sh in=<file> ref=<file,file,file...> pattern=out_%.fq outu=unmapped.fq
or
seal.sh in=<file> *.fasta.gz pattern=out_%.fq.gz outu=unmapped.fq.gz

Input may be fasta or fastq, compressed or uncompressed.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped 
fasta input, set in=stdin.fa.gz


Input parameters:
in=<file>           Main input. in=stdin.fq will pipe from stdin.
in2=<file>          Input for 2nd read of pairs in a different file.
ref=<file,file>     Comma-delimited list of reference files or directories.
                    Filenames may also be used without ref=, e.g. *.fa.
                    In addition to filenames, you may also use the keywords:
                    adapters, artifacts, phix, lambda, pjet, mtst, kapa.
literal=<seq,seq>   Comma-delimited list of literal reference sequences.
touppercase=f       (tuc) Change all bases upper-case.
interleaved=auto    (int) t/f overrides interleaved autodetection.
                    Must be set manually when streaming fastq input.
qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
reads=-1            If positive, quit after processing X reads or pairs.
copyundefined=f     (cu) Process non-AGCT IUPAC reference bases by making all
                    possible unambiguous copies.  Intended for short motifs
                    or adapter barcodes, as time/memory use is exponential.

Output parameters:
out=<file>          (outmatch) Write reads here that contain kmers matching
                    the reference. 'out=stdout.fq' will pipe to standard out.
out2=<file>         (outmatch2) Use this to write 2nd read of pairs to a 
                    different file.
outu=<file>         (outunmatched) Write reads here that do not contain kmers 
                    matching the database.
outu2=<file>        (outunmatched2) Use this to write 2nd read of pairs to a 
                    different file.
pattern=<file>      Use this to write reads to one stream per ref sequence
                    match, replacing the % character with the sequence name.
                    For example, pattern=%.fq for ref sequences named dog and 
                    cat would create dog.fq and cat.fq.
stats=<file>        Write statistics about which contaminants were detected.
refstats=<file>     Write statistics on a per-reference-file basis.
rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
dump=<file>         Dump kmer tables to a file, in fasta format.
nzo=t               Only write statistics about ref sequences with nonzero hits.
overwrite=t         (ow) Grant permission to overwrite files.
showspeed=t         (ss) 'f' suppresses display of processing speed.
ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
fastawrap=80        Length of lines in fasta output.
qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
statscolumns=5      (cols) Number of columns for stats output, 3 or 5.
                    5 includes base counts.
rename=f            Append matched reference sequence names to read headers.
addcount=t          If renaming, include the ref hit counts.
tophitonly=f        If renaming, only add the top hit.
refnames=f          Use names of reference files rather than scaffold IDs.
                    With multiple reference files, this is more efficient
                    than tracking statistics on a per-sequence basis.
trd=f               Truncate read and ref names at the first whitespace.
ordered=f           Set to true to output reads in same order as input.
kpt=t               (keepPairsTogether) Paired reads will always be assigned
                    to the same ref sequence.

Processing parameters:
k=31                Kmer length used for finding contaminants.  Contaminants 
                    shorter than k will not be found.  k must be at least 1.
rcomp=t             Look for reverse-complements of kmers in addition to 
                    forward kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to 
                    increase sensitivity in the presence of errors.  This may
                    also be set to a number, e.g. mm=3, to mask that many bp.
                    The default mm=t corresponds to mm=1 for odd-length kmers
                    and mm=2 for even-length kmers (as of v39.04), while
                    mm=f is always equivalent to mm=0.
minkmerhits=1       (mkh) A read needs at least this many kmer hits to be 
                    considered a match.
minkmerfraction=0.0 (mkf) A reads needs at least this fraction of its total
                    kmers to hit a ref, in order to be considered a match.
hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
                    Memory use is proportional to (3*K)^hdist.
qhdist=0            Hamming distance for query kmers; impacts speed, not memory.
editdistance=0      (edist) Maximum edit distance from ref kmers (subs and 
                    indels).  Memory use is proportional to (8*K)^edist.
forbidn=f           (fn) Forbids matching of read kmers containing N.  
                    By default, these will match a reference 'A' if hdist>0
                    or edist>0, to increase sensitivity.
match=all           Determines when to quit looking for kmer matches.  Values:
                         all:    Attempt to match all kmers in each read.
                         first:  Quit after the first matching kmer.
                         unique: Quit after the first uniquely matching kmer.
ambiguous=random    (ambig) Set behavior on ambiguously-mapped reads (with an
                    equal number of kmer matches to multiple sequences).
                         first:  Use the first best-matching sequence.
                         toss:   Consider unmapped.
                         random: Select one best-matching sequence randomly.
                         all:    Use all best-matching sequences.
genesets=f          Assign ambiguously-mapping reads to a newly created gene
                    set that they share for stats/rpkm output.  May be slow.
clearzone=0         (cz) Threshhold for ambiguity.  If the best match shares X 
                    kmers with the read, the read will be considered
                    also ambiguously mapped to any sequence sharing at least
                    [X minus clearzone] kmers.
czf=0.0             (clearzonefraction) If positive, the actual clearzone used
                    for a read with N total kmers will be max(cz, czf*N).
ecco=f              For overlapping paired reads only.  Performs error-
                    correction with BBMerge prior to kmer operations.

Containment parameters:
processcontainedref=f  Require a reference sequence to be fully contained by
                    an input sequence
storerefbases=f     Store reference bases so that ref containments can be
                    validated.  If this is set to false and processcontainedref
                    is true, then it will only require that the read share the
                    same number of bases as are present in the ref sequence.

Taxonomy parameters (only use when doing taxonomy):
tax=<file>          Output destination for taxonomy information.
taxtree=<file>      (tree) A serialized TaxTree (tree.taxtree.gz).
gi=<file>           A serialized GiTable (gitable.int1d.gz). Only needed if 
                    reference sequence names start with 'gi|'.
mincount=1          Only display taxa with at least this many hits.
maxnodes=-1         If positive, display at most this many top hits.
minlevel=subspecies Do not display nodes below this taxonomic level.
maxlevel=life       Do not display nodes above this taxonomic level.
Valid levels are subspecies, species, genus, family, order, class,
phylum, kingdom, domain, life

Speed and Memory parameters:
threads=auto        (t) Set number of threads to use; default is number of 
                    logical processors.
prealloc=f          Preallocate memory in table.  Allows faster table loading 
                    and more efficient memory usage, for a large reference.
monitor=f           Kill this process if CPU usage drops to zero for a long
                    time.  monitor=600,0.01 would kill after 600 seconds 
                    under 1% usage.
rskip=1             Skip reference kmers to reduce memory usage.
                    1 means use all, 2 means use every other kmer, etc.
qskip=1             Skip query kmers to increase speed.  1 means use all.
speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                    reads and reference.  Increases speed and reduces memory.
Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

Trimming/Masking parameters:
qtrim=f             Trim read ends to remove bases with quality below trimq.
                    Performed AFTER looking for kmers.  Values: 
                         t (trim both ends), 
                         f (neither end), 
                         r (right end only), 
                         l (left end only).
trimq=6             Regions with average quality BELOW this will be trimmed.
minlength=1         (ml) Reads shorter than this after trimming will be 
                    discarded.  Pairs will be discarded only if both are shorter.
maxlength=          Reads longer than this after trimming will be discarded.
                    Pairs will be discarded only if both are longer.
minavgquality=0     (maq) Reads with average quality (after trimming) below 
                    this will be discarded.
maqb=0              If positive, calculate maq from this many initial bases.
maxns=-1            If non-negative, reads with more Ns than this 
                    (after trimming) will be discarded.
forcetrimleft=0     (ftl) If positive, trim bases to the left of this position 
                    (exclusive, 0-based).
forcetrimright=0    (ftr) If positive, trim bases to the right of this position 
                    (exclusive, 0-based).
forcetrimright2=0   (ftr2) If positive, trim this many bases on the right end. 
forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
                    modulo this number.
restrictleft=0      If positive, only look for kmer matches in the 
                    leftmost X bases.
restrictright=0     If positive, only look for kmer matches in the 
                    rightmost X bases.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an 
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

seal() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP jgi.Seal $@"
	echo $CMD >&2
	eval $CMD
}

seal "$@"
