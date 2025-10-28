#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 20, 2024

Description:  RQCFilter2 is a revised version of RQCFilter that uses a common path for all dependencies.
The dependencies are available at http://portal.nersc.gov/dna/microbial/assembly/bushnell/RQCFilterData.tar

Performs quality-trimming, artifact removal, linker-trimming, adapter trimming, and spike-in removal using BBDuk.
Performs human/cat/dog/mouse/microbe removal using BBMap.
It requires 40 GB RAM for mousecatdoghuman, but only 1GB or so without them.

Usage:  rqcfilter2.sh in=<input file> path=<output directory> rqcfilterdata=<path to RQCFilterData directory>

Primary I/O parameters:
in=<file>           Input reads.
in2=<file>          Use this if 2nd read of pairs are in a different file.
path=null           Set to the directory to use for all output files.

Reference file path parameters:
rqcfilterdata=      Path to unzipped RQCFilterData directory.  Default is /global/projectb/sandbox/gaag/bbtools/RQCFilterData
ref=<file,file>     Comma-delimited list of additional reference files for filtering via BBDuk.

Output parameters:
scafstats=scaffoldStats.txt  Scaffold stats file name (how many reads matched which reference scaffold) .
kmerstats=kmerStats.txt      Kmer stats file name (duk-like output).
log=status.log               Progress log file name.
filelist=file-list.txt       List of output files.
stats=filterStats.txt        Overall stats file name.
stats2=filterStats2.txt      Better overall stats file name.
ihist=ihist_merge.txt        Insert size histogram name.  Set to null to skip merging.
outribo=ribo.fq.gz           Output for ribosomal reads, if removeribo=t.
reproduceName=reproduce.sh   Name of shellscript to reproduce these results.
usetmpdir=t                  Write temp files to TMPDIR.
tmpdir=                      Override TMPDIR.

Adapter trimming parameters:
trimhdist=1         Hamming distance used for trimming.
trimhdist2=         Hamming distance used for trimming with short kmers.  If unset, trimhdist will be used.
trimk=23            Kmer length for trimming stage.
mink=11             Minimum kmer length for short kmers when trimming.
trimfragadapter=t   Trim all known Illumina adapter sequences, including TruSeq and Nextera.
trimrnaadapter=f    Trim Illumina TruSeq-RNA adapters.
bisulfite=f         Currently, this trims the last 1bp from all reads after the adapter-trimming phase.
findadapters=t      For paired-end files, attempt to discover the adapter sequence with BBMerge and use that rather than a set of known adapters.
swift=f             Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.

Quality trimming parameters:
qtrim=f             Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers.
                    Values: rl (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10            Trim quality threshold.  Must also set qtrim for direction.
minlength=45        (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
mlf=0.333           (minlengthfraction) Reads shorter than this fraction of original length after trimming will be discarded.
minavgquality=5     (maq) Reads with average quality (before trimming) below this will be discarded.
maxns=0             Reads with more Ns than this will be discarded.
forcetrimmod=5      (ftm) If positive, right-trim length to be equal to zero, modulo this number.
forcetrimleft=-1    (ftl) If positive, trim bases to the left of this position
                    (exclusive, 0-based).
forcetrimright=-1   (ftr) If positive, trim bases to the right of this position
                    (exclusive, 0-based).
forcetrimright2=-1  (ftr2) If positive, trim this many bases on the right end.

Mapping parameters (for vertebrate contaminants):
mapk=14             Kmer length for mapping stage (9-15; longer is faster).
removehuman=f       (human) Remove human reads via mapping.
keephuman=f         Keep reads that map to human (or cat, dog, mouse) rather than removing them.
removedog=f         (dog) Remove dog reads via mapping.
removecat=f         (cat) Remove cat reads via mapping.
removemouse=f       (mouse) Remove mouse reads via mapping.
aggressivehuman=f   Aggressively remove human reads (and cat/dog/mouse) using unmasked references.
aggressivemicrobe=f Aggressively microbial contaminant reads using unmasked references.
aggressive=f        Set both aggressivehuman and aggressivemicrobe at once.
mapref=             Remove contaminants by mapping to this fasta file (or comma-delimited list).

Bloom filter parameters (for vertebrate mapping):
bloom=t             Use a Bloom filter to accelerate mapping.
bloomminreads=4m   Disable Bloom filter if there are fewer than this many reads.
bloomk=29           Kmer length for Bloom filter
bloomhashes=1       Number of hashes for the Bloom filter.
bloomminhits=6      Minimum consecutive hits to consider a read as matching.
bloomserial=t       Use the serialized Bloom filter for greater loading speed.
                    This will use the default Bloom filter parameters.

Microbial contaminant removal parameters:
detectmicrobes=f    Detect common microbes, but don't remove them.  Use this OR removemicrobes, not both.
removemicrobes=f    (microbes) Remove common contaminant microbial reads via mapping, and place them in a separate file.
taxlist=            (tax) Remove these taxa from the database before filtering.  Typically, this would be the organism name or NCBI ID, or a comma-delimited list.  Organism names should have underscores instead of spaces, such as Escherichia_coli.
taxlevel=order      (level) Level to remove.  For example, 'phylum' would remove everything in the same phylum as entries in the taxlist.
taxtree=auto        (tree) Override location of the TaxTree file.
gitable=auto        Override location of the gitable file.
loadgitable=f       Controls whether gi numbers may be used for taxonomy.
microberef=         Path to fasta file of microbes.
microbebuild=1      Chooses which masking was used.  1 is most stringent and should be used for bacteria.  Eukaryotes should use 3.

Extended microbial contaminant parameters:
detectmicrobes2=f   (detectothermicrobes) Detect an extended set of microbes that are currently being screened.  This can be used in conjunction with removemicrobes.

Filtering parameters (for artificial and genomic contaminants):
skipfilter=f        Skip this phase.  Not recommended.
filterpolya=f       Remove reads containing poly-A sequence (for RNA-seq).
filterpolyg=0       Remove reads that start with a G polymer at least this long (0 disables).
trimpolyg=6         Trim reads that start or end with a G polymer at least this long (0 disables).
phix=t              Remove reads containing phiX kmers.
lambda=f            Remove reads containing Lambda phage kmers.
pjet=t              Remove reads containing PJET kmers.
sip=f               Remove SIP spikeins.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors.
maxbadkmers=0       (mbk) Reads with more than this many contaminant kmers will be discarded.
filterhdist=1       Hamming distance used for filtering.
filterqhdist=1      Query hamming distance used for filtering.
copyundefined=f     (cu) Match all possible bases for sequences containing degerate IUPAC symbols.
entropy=f           Remove low-complexity reads.  The threshold can be specified by e.g entropy=0.4; default is 0.42 if enabled.
entropyk=2          Kmer length to use for entropy calculation.
entropywindow=40    Window size to use for entropy calculation.

Spikein removal/quantification parameters:
mtst=f              Remove mtst.
kapa=t              Remove and quantify kapa.
spikeink=31         Kmer length for spikein removal.
spikeinhdist=0      Hamming distance for spikein removal.
spikeinref=         Additional references for spikein removal (comma-delimited list).

Ribosomal filtering parameters:
ribohdist=1         Hamming distance used for rRNA removal.
riboedist=0         Edit distance used for rRNA removal.
removeribo=f        (ribo) Remove ribosomal reads via kmer-matching, and place them in a separate file.

Organelle filtering parameters:
chloromap=f         Remove chloroplast reads by mapping to this organism's chloroplast.
mitomap=f           Remove mitochondrial reads by mapping to this organism's mitochondria.
ribomap=f           Remove ribosomal reads by mapping to this organism's ribosomes.
NOTE: organism TaxID should be specified in taxlist, and taxlevel should be set to genus or species.

FilterByTile parameters:
filterbytile=t      Run FilterByTile to remove reads from low-quality parts of the flowcell.
tiledump=           Set this to the tiledump of the full lane (recommended).

Recalibration parameters:
recalibrate=t       Recalibrate quality scores based on PhiX alignment.
phixsam=            Set this to the aligned PhiX data for the lane (required).
quantize=2          Quantize the quality scores to reduce file size, using this divisor.
                    2 reduces size by roughly 25%.  Disabled if recalibrate=f.
                    Quantization happens AFTER all the quality-related steps.

Polyfilter parameters:
polyfilter=GC       Remove reads with homopolymers of these subunits.
                    Set polyfilter=null to disable.

Clumpify parameters:
clumpify=f          Run clumpify; all deduplication flags require this.
dedupe=f            Remove duplicate reads; all deduplication flags require this.
opticaldupes=f      Remove optical duplicates (Clumpify optical flag).
edgedupes=f         Remove tile-edge duplicates (Clumpify spany and adjacent flags).
dpasses=1           Use this many deduplication passes.
dsubs=2             Allow this many substitutions between duplicates.
ddist=40            Remove optical/edge duplicates within this distance.
lowcomplexity=f     Set to true for low-complexity libraries such as RNA-seq to improve estimation of memory requirements.
clumpifytmpdir=f    Use TMPDIR for clumpify temp files.
clumpifygroups=-1   If positive, force Clumpify to use this many groups.
*** For NextSeq, the recommended deduplication flags are: clumpify dedupe edgedupes
*** For NovaSeq, the recommended deduplication flags are: clumpify dedupe opticaldupes ddist=12000
*** For HiSeq, the recommended deduplication flags are: clumpify dedupe opticaldupes

Sketch parameters:
sketch=t            Run SendSketch on 2M read pairs.
silvalocal=t        Use the local flag for Silva (requires running RQCFilter on NERSC).
sketchreads=1m      Number of read pairs to sketch.
sketchsamplerate=1  Samplerate for SendSketch.
sketchminprob=0.2   Minprob for SendSketch.
sketchdb=nt,refseq,silva  Servers to use for SendSketch.

Other processing parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
library=frag        Set to 'frag', 'clip', 'lfpe', or 'clrs'.
filterk=31          Kmer length for filtering stage.
rcomp=t             Look for reverse-complements of kmers in addition to forward kmers.
nexteralmp=f        Split into different files based on Nextera LMP junction sequence.  Only for Nextera LMP, not normal Nextera.
extend=f            Extend reads during merging to allow insert size estimation of non-overlapping reads.
pigz=t              Use pigz for compression.
unpigz=t            Use pigz for decompression.
khist=f             Set to true to generate a kmer-frequency histogram of the output data.
merge=t             Set to false to skip generation of insert size histogram.

Header-specific parameters:  (NOTE - Be sure to disable these if the reads have improper headers, like SRA data.)
chastityfilter=t    Remove reads failing chastity filter.
barcodefilter=f     Crash when improper barcodes are discovered.  Set to 'f' to disable,
                    't' to remove improper barcodes, or 'crash' to crash if they are discovered.
barcodes=           A comma-delimited list of barcodes or files of barcodes.
filterbytile        Also needs to be disabled for SRA data.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

*****   All additional parameters supported by BBDuk may also be used, and will be passed directly to BBDuk   *****

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

z="-Xmx40g"
z2="-Xms40g"
set=0
export TZ="America/Los_Angeles"

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
	freeRam 39200m 84

	if [[ $NSLOTS == 8 ]]; then
		RAM=39200
	fi
	
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}

calcXmx "$@"


rqcfilter() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
		export TZ="America/Los_Angeles" 
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load pigz
		export TZ="America/Los_Angeles" 
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load pigz
	fi
	local CMD="java $EA $SIMD $EOOM $z $z2 $JNI -cp $CP jgi.RQCFilter2 jni=t $@"
	echo $CMD >&2
	eval $CMD
}

rqcfilter "$@"
