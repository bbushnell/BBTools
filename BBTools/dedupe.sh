#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Jonathan Rood
Last modified February 19, 2020

Description:  Accepts one or more files containing sets of sequences (reads or scaffolds).
Removes duplicate sequences, which may be specified to be exact matches, subsequences, or sequences within some percent identity.
Can also find overlapping sequences and group them into clusters.
Please read bbmap/docs/guides/DedupeGuide.txt for more information.

Usage:     dedupe.sh in=<file or stdin> out=<file or stdout>

An example of running Dedupe for clustering short reads:
dedupe.sh in=x.fq am=f ac=f fo c pc rnc=f mcs=4 mo=100 s=1 pto cc qin=33 csf=stats.txt pattern=cluster_%.fq dot=graph.dot

Input may be fasta or fastq, compressed or uncompressed.
Output may be stdout or a file.  With no output parameter, data will be written to stdout.
If 'out=null', there will be no output, but statistics will still be printed.
You can also use 'dedupe <infile> <outfile>' without the 'in=' and 'out='.

I/O parameters:
in=<file,file>        A single file or a comma-delimited list of files.
out=<file>            Destination for all output contigs.
pattern=<file>        Clusters will be written to individual files, where the '%' symbol in the pattern is replaced by cluster number.
outd=<file>           Optional; removed duplicates will go here.
csf=<file>            (clusterstatsfile) Write a list of cluster names and sizes.
dot=<file>            (graph) Write a graph in dot format.  Requires 'fo' and 'pc' flags.
threads=auto          (t) Set number of threads to use; default is number of logical processors.
overwrite=t           (ow) Set to false to force the program to abort rather than overwrite an existing file.
showspeed=t           (ss) Set to 'f' to suppress display of processing speed.
minscaf=0             (ms) Ignore contigs/scaffolds shorter than this.
interleaved=auto      If true, forces fastq input to be paired and interleaved.
ziplevel=2            Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.

Output format parameters:
storename=t           (sn) Store scaffold names (set false to save memory).
#addpairnum=f         Add .1 and .2 to numeric id of read1 and read2.
storequality=t        (sq) Store quality values for fastq assemblies (set false to save memory).
uniquenames=t         (un) Ensure all output scaffolds have unique names.  Uses more memory.
mergenames=f          When a sequence absorbs another, concatenate their headers.
mergedelimiter=>      Delimiter between merged headers.  Can be a symbol name like greaterthan.
numbergraphnodes=t    (ngn) Label dot graph nodes with read numbers rather than read names.
sort=f                Sort output (otherwise it will be random).  Options:
                         length:  Sort by length
                         quality: Sort by quality
                         name:    Sort by name
                         id:      Sort by input order
ascending=f           Sort in ascending order.
ordered=f             Output sequences in input order.  Equivalent to sort=id ascending.
renameclusters=f      (rnc) Rename contigs to indicate which cluster they are in.
printlengthinedges=f  (ple) Print the length of contigs in edges.

Processing parameters:
absorbrc=t            (arc) Absorb reverse-complements as well as normal orientation.
absorbmatch=t         (am) Absorb exact matches of contigs.
absorbcontainment=t   (ac) Absorb full containments of contigs.
#absorboverlap=f      (ao) Absorb (merge) non-contained overlaps of contigs (TODO).
findoverlap=f         (fo) Find overlaps between contigs (containments and non-containments).  Necessary for clustering.
uniqueonly=f          (uo) If true, all copies of duplicate reads will be discarded, rather than keeping 1.
rmn=f                 (requirematchingnames) If true, both names and sequence must match.
usejni=f              (jni) Do alignments in C code, which is faster, if an edit distance is allowed.
                      This will require compiling the C code; details are in /jni/README.txt.

Subset parameters:
subsetcount=1         (sstc) Number of subsets used to process the data; higher uses less memory.
subset=0              (sst) Only process reads whose ((ID%subsetcount)==subset).

Clustering parameters:
cluster=f             (c) Group overlapping contigs into clusters.
pto=f                 (preventtransitiveoverlaps) Do not look for new edges between nodes in the same cluster.
minclustersize=1      (mcs) Do not output clusters smaller than this.
pbr=f                 (pickbestrepresentative) Only output the single highest-quality read per cluster.

Cluster postprocessing parameters:
processclusters=f     (pc) Run the cluster processing phase, which performs the selected operations in this category.
                      For example, pc AND cc must be enabled to perform cc.
fixmultijoins=t       (fmj) Remove redundant overlaps between the same two contigs.
removecycles=t        (rc) Remove all cycles so clusters form trees.
cc=t                  (canonicizeclusters) Flip contigs so clusters have a single orientation.
fcc=f                 (fixcanoncontradictions) Truncate graph at nodes with canonization disputes.
foc=f                 (fixoffsetcontradictions) Truncate graph at nodes with offset disputes.
mst=f                 (maxspanningtree) Remove cyclic edges, leaving only the longest edges that form a tree.

Overlap Detection Parameters
exact=t               (ex) Only allow exact symbol matches.  When false, an 'N' will match any symbol.
touppercase=t         (tuc) Convert input bases to upper-case; otherwise, lower-case will not match.
maxsubs=0             (s) Allow up to this many mismatches (substitutions only, no indels).  May be set higher than maxedits.
maxedits=0            (e) Allow up to this many edits (subs or indels).  Higher is slower.
minidentity=100       (mid) Absorb contained sequences with percent identity of at least this (includes indels).
minlengthpercent=0    (mlp) Smaller contig must be at least this percent of larger contig's length to be absorbed.
minoverlappercent=0   (mop) Overlap must be at least this percent of smaller contig's length to cluster and merge.
minoverlap=200        (mo) Overlap must be at least this long to cluster and merge.
depthratio=0          (dr) When non-zero, overlaps will only be formed between reads with a depth ratio of at most this.
                      Should be above 1.  Depth is determined by parsing the read names; this information can be added
                      by running KmerNormalize (khist.sh, bbnorm.sh, or ecc.sh) with the flag 'rename'
k=31                  Seed length used for finding containments and overlaps.  Anything shorter than k will not be found.
numaffixmaps=1        (nam) Number of prefixes/suffixes to index per contig. Higher is more sensitive, if edits are allowed.
hashns=f              Set to true to search for matches using kmers containing Ns.  Can lead to extreme slowdown in some cases.
#ignoreaffix1=f       (ia1) Ignore first affix (for testing).
#storesuffix=f        (ss) Store suffix as well as prefix.  Automatically set to true when doing inexact matches.

Other Parameters
qtrim=f               Set to qtrim=rl to trim leading and trailing Ns.
trimq=6               Quality trim level.
forcetrimleft=-1      (ftl) If positive, trim bases to the left of this position (exclusive, 0-based).
forcetrimright=-1     (ftr) If positive, trim bases to the right of this position (exclusive, 0-based).

Note on Proteins / Amino Acids
Dedupe supports amino acid space via the 'amino' flag.  This also changes the default kmer length to 10.
In amino acid mode, all flags related to canonicity and reverse-complementation are disabled,
and nam (numaffixmaps) is currently limited to 2 per tip.

Java Parameters:
-Xmx                  This will set Java's memory usage, overriding autodetection.
                      -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom                 This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da                   Disable assertions.

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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

dedupe() {
	local CMD="java $JNI $EA $EOOM $z $z2 -cp $CP jgi.Dedupe $@"
	echo $CMD >&2
	eval $CMD
}

dedupe "$@"
