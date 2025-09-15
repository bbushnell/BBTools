This is a readme for BBMap.  However, it has not been maintained and is superceded by the information in the shellscript, bbmap.sh.

Basic Syntax:

(Using shellscript, under Unix, which autodetects RAM to set -Xmx parameter.  You can also include a flag like '-Xmx31g' in the shellscript arguments to set RAM usage.)
To index:
bbmap.sh ref=<reference.fa>
To map:
bbmap.sh in=<reads.fq> out=<mapped.sam>

(without shellscript)
To index:
java -ea -Xmx31g -cp <PATH> align2.BBMap ref=<reference.fa>
To map:
java -ea -Xmx31g -cp <PATH> align2.BBMap in=<reads.fq> out=<mapped.sam>

...where "<PATH>" should indicate the path to the directory containing all the source code directories; e.g. "/usr/bin/bbmap/current"

Please note, the reference is only needed for building the index the first time; subsequently, just specify the build number which corresponds to that reference.
So for example the first time you map to e.coli you might specify "ref=ecoli_reference.fa build=3"; after that, just specify "build=3".
The index files would then be stored in ./ref/genome/3/ and ./ref/index/3/
Also, the -Xmx parameter should specify approximately 85% of the physical memory of the target machine; so, 21G for a 24GB node.  The process needs approximately 8 bytes per reference base (plus a several hundred MB overhead).


Advanced Syntax:


Indexing Parameters (required when building the index):
path=<.>        	Base directory to store index files.  Default is the local directory.  The index will always be placed in a subdirectory "ref".
ref=<ref.fasta> 	Use this file to build the index.  Needs to be specified only once; subsequently, the build number should be used.
build=<1>		Write the index to this location (build=1 would be stored in /ref/genome/1/ and /ref/index/1/).  Can be any integer.  This parameter defaults to 1, but using additional numbers allows multiple references to be indexed in the same directory.
k=<13>          	Use length 13 kmers for indexing.  Suggested values are 9-15, with lower typically being slower and more accurate.  13 is usually optimal.  14 is better for RNA-SEQ and very large references >4GB; 12 is better for PacBio and cross-species mapping.
midpad=<300>		Put this many "N" in between scaffolds when making the index.  300 is fine for metagenomes with millions of contigs; for a finished genome like human with 25 scaffolds, this should be set to 100000+ to prevent cross-scaffold mapping.
startpad=<8000> 	Put this many "N" at the beginning of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
stoppad=<8000>		Put this many "N" at the end of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
minscaf=<1>		Do not include scaffolds shorter than this when generating index.  Useful for assemblies with millions of fairly worthless unscaffolded contigs under 100bp.  There's no reason to make this shorter than the kmer length.
usemodulo=<f>		Throw away ~80% of kmers based on their remainder modulo a number.  Reduces memory usage by around 50%, and reduces sensitivity slightly.  Must be specified when indexing and when mapping.


Input Parameters:
path=<.>		Base directory to read index files.
build=<1>		Use the index at this location (same as when indexing).
in=<reads.fq>		Use this as the input file for reads.  Also accepts fasta.  "in=sequential length=200" will break a genome into 200bp pieces and map them to itself.  "in=stdin" will accept piped input.  The format of piped input can be specified with e.g. "in=stdin.fq.gz" or "in=stdin.fa"; default is uncompressed fastq.
in2=<reads2.fq> 	Run mapping paired, with reads2 in the file "reads2.fq"
			NOTE:  As a shorthand, "in=reads#.fq" is equivalent to "in=reads1.fq in2=reads2.fq"
interleaved=<auto>	Or "int". Set to "true" to run mapping paired, forcing the reads to be considered interleaved from a single input file.  By default the reader will try to determine whether a file is interleaved based on the read names; so if you don't want this, set interleaved=false.
qin=<auto>       	Set to 33 or 64 to specify input quality value ASCII offset.
fastareadlen=<500>	If fasta is used for input, breaks the fasta file up into reads of about this length.  Useful if you want to map one reference against another, since BBMap currently has internal buffers limited to 500bp.  I can change this easily if desired.
fastaminread=<1>	Ignore fasta reads shorter than this.  Useful if, say, you set fastareadlen=500, and get a length 518 read; this will be broken into a 500bp read and an 18bp read.  But it's not usually worth mapping the 18bp read, which will often be ambiguous.
maxlen=<0>        	Break long fastq reads into pieces of this length.
minlen=<0>       	Throw away remainder of read that is shorter than this.
fakequality=<-1>	Set to a positive number 1-50 to generate fake quality strings for fasta input reads.  Less than one turns this function off.
blacklist=<a.fa,b.fa>	Set a list of comma-delimited fasta files.  Any read mapped to a scaffold name in these files will be considered "blacklisted" and can be handled differently by using the "outm", "outb", and "outputblacklisted" flags.  The blacklist fasta files should also be merged with other fasta files to make a single combined fasta file; this combined file should be specified with the "ref=" flag when indexing.
touppercase=<f>		Set true to convert lowercase read bases to upper case.  This is required if any reads have lowercase letters (which real reads should never have).


Sampling Parameters:
reads=<-1>		Process at most N reads, then stop.  Useful for benchmarking.  A negative number will use all reads.
samplerate=<1.0>	Set to a fraction of 1 if you want to randomly sample reads.  For example, samplerate=0.25 would randomly use a quarter of the reads and ignore the rest.  Useful for huge datasets where all you want to know is the % mapped.
sampleseed=<1>		Set to the RNG seed for random sampling.  If this is set to a negative number, a random seed is used; for positive numbers, the number itself is the seed.  Since the default is 1, this is deterministic unless you explicitly change it to a negative number.	
idmodulo=<1>		Set to a higher number if you want to map only every Nth read (for sampling huge datasets).


Mapping Parameters:
fast=<f>		The fast flag is a macro.  It will set many other paramters so that BBMap will run much faster, at slightly reduced sensitivity for most applications.  Not recommended for RNAseq, cross-species alignment, or other situations where long deletions or low identity matches are expected.
minratio=<0.56>		Alignment sensitivity as a fraction of a read's max possible mapping score.  Lower is slower and more sensitive but gives more false positives.  Ranges from 0 (very bad alignment) to 1 (perfect alignment only).  Default varies between BBMap versions. 
minidentity=<>		Or "minid".  Use this flag to set minratio more easily.  If you set minid=0.9, for example, minratio will be set to a value that will be APPROXIMATELY equivalent to 90% identity alignments.
minapproxhits=<1>	Controls minimum number of seed hits to examine a site.  Higher is less accurate but faster (on large genomes).  2 is maybe 2.5x as fast and 3 is maybe 5x as fast on a genome with several gigabases.  Does not speed up genomes under 100MB or so very much.
padding=<4>		Sets extra padding for slow-aligning.  Higher numbers are more accurate for indels near the tips of reads, but slower.
tipsearch=<100>		Controls how far to look for possible deletions near tips of reads by brute force.  tipsearch=0 disables this function.  Higher is more accurate.
maxindel=<16000>	Sets the maximum size of indels allowed during the quick mapping phase.  Set higher (~100,000) for RNA-SEQ and lower (~20) for large assemblies with mostly very short contigs.  Lower is faster.
strictmaxindel=<f>	Set to true to disallow mappings with indels longer than maxindel.  Alternately, for an integer X, 'strictmaxindel=X' is equivalent to the pair of flags 'strictmaxindel=t maxindel=X'.
pairlen=<32000>  	Maximum distance between mates allowed for pairing.
requirecorrectstrand=<t>	Or "rcs".  Requires correct strand orientation when pairing reads.  Please set this to false for long mate pair libraries!
samestrandpairs=<f>	Or "ssp".  Defines correct strand orientation when pairing reads.  Default is false, meaning opposite strands, as in Illumina fragment libraries.  "ssp=true" mode is not fully tested.
killbadpairs=<f>	Or "kbp".  When true, if a read pair is mapped with an inappropriate insert size or orientation, the read with the lower mapping quality is marked unmapped.
rcompmate=<f>		***TODO*** Set to true if you wish the mate of paired reads to be reverse-complemented prior to mapping (to allow better pairing of same-strand pair libraries).
kfilter=<-1>		If set to a positive number X, all potential mapping locatiosn that do not have X contiguous perfect matches with the read will be ignored.  So, reads that map with "kfilter=51" are assured to have at least 51 contiguous bases that match the reference.  Useful for mapping to assemblies generated by a De Bruijn graph assembly that used a kmer length of X, so that you know which reads were actually used in the assembly.
threads=<?>		Or "t".  Set number of threads.  Default is # of logical cores.  The total number of active threads will be higher than this, because input and output are in seperate threads.
perfectmode=<f>		Only accept perfect mappings.  Everything goes much faster.  
semiperfectmode=<f>	Only accept perfect or "semiperfect" mappings.  Semiperfect means there are no mismatches of defined bases, but up to half of the reference is 'N' (to allow mapping to the edge of a contig).
rescue=<t>		Controls whether paired may be rescued by searching near the mapping location of a mate.  Increases accuracy, with usually a minor speed penalty.
expectedsites=<1>	For BBMapPacBioSkimmer only, sets the expected number of correct mapping sites in the target reference.  Useful if you are mapping reads to other reads with some known coverage depth.
msa=<>			Advanced option, not recommended.  Set classname of MSA to use.
bandwidth=0		Or "bw".  When above zero, restricts alignment band to this width.  Runs faster, but with reduced accuracy for reads with many or long indels.
bandwidthratio=0	Or "bwr".  When above zero, restricts alignment band to this fraction of a read's length.  Runs faster, but with reduced accuracy for reads with many or long indels.
usequality=<t>		Or "uq".  Set to false to ignore quality values when mapping.  This will allow very low quality reads to be attempted to be mapped rather than discarded.
keepbadkeys=<f>		Or "kbk".  With kbk=false (default), read keys (kmers) have their probability of being incorrect evaluated from quality scores, and keys with a 94%+ chance of being wrong are discarded.  This increases both speed and accuracy.
usejni=<f>		Or "jni".  Do alignments in C code, which is faster.  Requires first compiling the C code; details are in /jni/README.txt.  This will produce identical output.
maxsites2=<800>		Don't analyze (or print) more than this many alignments per read.
minaveragequality=<0>	(maq) Discard reads with average quality below this.

Post-Filtering Parameters:

idfilter=0              Different than "minid".  No alignments will be allowed with an identity score lower than this value.  This filter occurs at the very end and is unrelated to minratio, and has no effect on speed unless set to 1.  Range is 0-1.
subfilter=-1            Ban alignments with more than this many substitutions.
insfilter=-1            Ban alignments with more than this many insertions.
delfilter=-1            Ban alignments with more than this many deletions.
indelfilter=-1          Ban alignments with more than this many indels.
editfilter=-1           Ban alignments with more than this many edits.
inslenfilter=-1         Ban alignments with an insertion longer than this.
dellenfilter=-1         Ban alignments with a deletion longer than this.

Output Parameters:
out=<outfile.sam>	Write output to this file.  If out=null, output is suppressed.  If you want to output paired reads to paired files, use a "#" symbol, like out=mapped#.sam.  Then reads1 will go to mapped1.sam and reads2 will go to mapped2.sam. (NOTE: split output currently diabled for .sam format, but allowed for native .txt format).  To print to standard out, use "out=stdout"
outm=<>			Write only mapped reads to this file (excluding blacklisted reads, if any).
outu=<>			Write only unmapped reads to this file.
outb=<>			Write only blacklisted reads to this file.  If a pair has one end mapped to a non-blacklisted scaffold, it will NOT go to this file. (see: blacklist)
out2=<>			If you set out2, outu2, outm2, or outb2, the second read in each pair will go to this file.  Not currently allowed for SAM format, but OK for others (such as fasta, fastq, bread).
overwrite=<f>		Or "ow".  Overwrite output file if it exists, instead of aborting.
append=<f>		Or "app".  Append to output file if it exists, instead of aborting.
ambiguous=<best>	Or "ambig". Sets how to handle ambiguous reads.  "first" or "best" uses the first encountered best site (fastest).  "all" returns all best sites.  "random" selects a random site from all of the best sites (does not yet work with paired-ends).  "toss" discards all sites and considers the read unmapped (same as discardambiguous=true).  Note that for all options (aside from toss) ambiguous reads in SAM format will have the extra field "XT:A:R" while unambiguous reads will have "XT:A:U".
ambiguous2=<best>	(for BBSplit only) Or "ambig2". Only for splitter mode.  Ambiguous2 strictly refers to any read that maps to more than one reference set, regardless of whether it has multiple mappings within a reference set.  This may be set to "best" (aka "first"), in which case the read will be written only to the first reference to which it has a best mapping; "all", in which case a read will be written to outputs for all references to which it maps; "toss", in which case it will be considered unmapped; or "split", in which case it will be written to a special output file with the prefix "AMBIGUOUS_" (one per reference).
outputunmapped=<t>	Outputs unmapped reads to primary output stream (otherwise they are dropped).
outputblacklisted=<t>	Outputs blacklisted reads to primary output stream (otherwise they are dropped).
ordered=<f>		Set to true if you want reads to be output in the same order they were input.  This takes more memory, and can be slower, due to buffering in multithreaded execution.  Not needed for singlethreaded execution.
ziplevel=<2>		Sets output compression level, from 1 (fast) to 9 (slow).  I/O is multithreaded, and thus faster when writing paired reads to two files rather than one interleaved file.
nodisk=<f>		"true" will not write the index to disk, and may load slightly faster.   Prevents collisions between multiple bbmap instances writing indexes to the same location at the same time.
usegzip=<f>		If gzip is installed, output file compression is done with a gzip subprocess instead of with Java's native deflate method.  Can be faster when set to true.  The output file must end in a compressed file extension for this to have effect.
usegunzip=<f>		If gzip is installed, input file decompression is done with a gzip subprocess instead of with Java's native inflate method.  Can be faster when set to true.
pigz=<f>          	Spawn a pigz (parallel gzip) process for faster compression than Java or gzip.  Requires pigz to be installed.
unpigz=<f>        	Spawn a pigz process for faster decompression than Java or gzip.  Requires pigz to be installed.
bamscript=<filename>	(bs for short) Writes a shell script to <filename> with the command line to translate the sam output of BBMap into a sorted bam file, assuming you have samtools in your path.
maxsites=<5>		Sets maximum alignments to print per read, if secondary alignments are allowed.  Currently secondary alignments may lack cigar strings.
secondary=<f>		Print secondary alignments.
sssr=<0.95>  		(secondarysitescoreratio) Print only secondary alignments with score of at least this fraction of primary.
ssao=<f>     		(secondarysiteasambiguousonly) Only print secondary alignments for ambiguously-mapped reads.
quickmatch=<f>		Generate cigar strings during the initial alignment (before the best site is known).  Currently, this must be enabled to generate cigar strings for secondary alignments.  It increases overall speed but may in some very rare cases yield inferior alignments due to less padding.
local=<f>          	Output local alignments instead of global alignments.  The mapping will still be based on the best global alignment, but the mapping score, cigar string, and mapping coordinate will reflect a local alignment (using the same affine matrix as the global alignment).
sortscaffolds=<f> 	Sort scaffolds alphabetically in SAM headers to allow easier comparisons with Tophat (in cuffdif, etc).  Default is in same order as source fasta.
trimreaddescriptions=<f>	(trd) Truncate read names at the first whitespace, assuming that the remaineder is a comment or description.
machineout=<f>    	Set to true to output statistics in machine-friendly 'key=value' format.
forcesectionname=<f>	All fasta reads get an _# at the end of their name.  The number is 1 for the first shred and continues ascending.


Sam settings and flags:
samversion=<1.4>  	SAM specification version. Set to 1.3 for cigar strings with 'M' or 1.4 for cigar strings with '=' and 'X'.  Samtools 0.1.18 and earlier are incompatible with sam format version 1.4 and greater.
saa=<t>           	(secondaryalignmentasterisks) Use asterisks instead of bases for sam secondary alignments.
cigar=<t>		Generate cigar strings (for bread format, this means match strings).  cigar=false is faster.  "cigar=" is synonymous with "match=".  This must be enabled if match/insertion/deletion/substitution statistics are desired, but the program will run faster with cigar strings disabled.
keepnames=<f>		Retain original names of paired reads, rather than ensuring both reads have the same name when written in sam format by renaming read2 to the same as read1.  If this is set to true then the output may not be sam compliant.
mdtag=<f>		Generate MD tags for SAM files.  Requires that cigar=true.  I do not recommend generating MD tags for RNASEQ or other data where long deletions are expected because they will be incredibly long.
xstag=<f>		Generate XS (strand) tags for Cufflinks.  This should be used with a stranded RNA-seq protocol.
xmtag=<t>		Generate XM tag.  Indicates number of best alignments.  May only work correctly with ambig=all.
nhtag=<f>             	Write NH tags.
intronlen=<999999999>	Set to a lower number like 10 to change 'D' to 'N' in cigar strings for deletions of at least that length.  This is used by Cufflinks; 'N' implies an intron while 'D' implies a deletion, but they are otherwise identical.
stoptag=<f>		Allows generation of custom SAM tag YS:i:<read stop location>
idtag=<f>		Allows generation of custom SAM tag YI:f:<percent identity>
scoretag=<f>		Allows generation of custom SAM tag YR:i:<raw mapping score>
inserttag=<f>		Write a tag indicating insert size, prefixed by X8:Z:
rgid=<>     		Set readgroup ID.  All other readgroup fields can be set similarly, with the flag rgXX=value.
noheader=<f>		Suppress generation of output header lines.


Statistics and Histogram Parameters:
showprogress=<f>	Set to true to print out a '.' once per million reads processed.  You can also change the interval with e.g. showprogress=20000.
qhist=<file>		Output a per-base average quality histogram to <file>.
aqhist=<file>		Write histogram of average read quality to <file>.
bqhist=<file>		Write a quality histogram designed for box plots to <file>.
obqhist=<file>		Write histogram of overall base counts per quality score to <file>.
qahist=<file>		Quality accuracy histogram; correlates claimed phred quality score with observed quality based on substitution, insertion, and deletion rates.
mhist=<file>		Output a per-base match histogram to <file>.  Requires cigar strings to be enabled.  The columns give fraction of bases at each position having each match string operation: match, substitution, deletion, insertion, N, or other.
ihist=<file>		Output a per-read-pair insert size histogram to <file>.
bhist=<file>		Output a per-base composition histogram to <file>.
indelhist=<file>  	Output an indel length histogram.
lhist=<file>      	Output a read length histogram.
ehist=<file>      	Output an errors-per-read histogram.
gchist=<file>     	Output a gc content histogram.
gchistbins=<100>	(gcbins) Set the number of bins in the gc content histogram.
idhist=<file>    	Write a percent identity histogram.
idhistbins=<100>	(idbins) Set the number of bins in the identity histogram.
scafstats=<file>	Track mapping statistics per scaffold, and output to <file>.
refstats=<file>		For BBSplitter, enable or disable tracking of read mapping statistics on a per-reference-set basis, and output to <file>.
verbosestats=<0>	From 0-3; higher numbers will print more information about internal program counters.
printunmappedcount=<f>	Set true to print the count of reads that were unmapped.  For paired reads this only includes reads whose mate was also unmapped.


Coverage output parameters (these may reduce speed and use more RAM):
covstats=<file>		Per-scaffold coverage info.
covhist=<file>		Histogram of # occurrences of each depth level.
basecov=<file>		Coverage per base location.
bincov=<file>		Print binned coverage per location (one line per X bases).
covbinsize=1000		Set the binsize for binned coverage output.
nzo=f			Only print scaffolds with nonzero coverage.
twocolumn=f		Change to true to print only ID and Avg_fold instead of all 6 columns to the 'out=' file.
32bit=f			Set to true if you need per-base coverage over 64k.
bitset=f		Store coverage data in BitSets.
arrays=t		Store coverage data in Arrays.
ksb=t			Keep residual bins shorter than binsize.
strandedcov=f		Track coverage for plus and minus strand independently.  Requires a # symbol in coverage output filenames which will be replaced by 1 for plus strand and 2 for minus strand.
startcov=f		Only track start positions of reads.
concisecov=f		Write basecov in a more concise format.


Trimming Parameters:
qtrim=<f>		Options are false, left, right, or both.  Allows quality-trimming of read ends before mapping.
			false: Disable trimming.
			left (l): Trim left (leading) end only.
			right (r): Trim right (trailing) end only.  This is the end with lower quality many platforms.
			both (lr): Trim both ends.
trimq=<5>		Set the quality cutoff.  Bases will be trimmed until there are 2 consecutive bases with quality GREATER than this value; default is 5.  If the read is from fasta and has no quality socres, Ns will be trimmed instead, as long as this is set to at least 1.
untrim=<f>		Untrim the read after mapping, restoring the trimmed bases.  The mapping position will be corrected (if necessary) and the restored bases will be classified as soft-clipped in the cigar string.


Java Parameters:
-Xmx       		If running from the shellscript, include it with the rest of the arguments and it will be passed to Java to set memory usage, overriding the shellscript's automatic memory detection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max allowed is typically 85% of physical memory.
-da			Disable assertions.  Alternative is -ea which is the default.


Splitting Parameters:
The splitter is invoked by calling bbsplit.sh (or align2.BBSplitter) instead of bbmap.sh, for the indexing phase.  It allows combining multiple references and outputting reads to different files depending on which one they mapped to best.  The order in which references are specified is important in cases of ambiguous mappings; when a read has 2 identically-scoring mapping locations from different references, it will be mapped to the first reference.
All parameters are the same as BBMap with the exception of the ones listed below.  You can still use "outu=" to capture unmapped reads.
ref_<name>=<fasta files>	Defines a named set of organisms with a single fasta file or list.  For example, ref_a=foo.fa,bar.fa defines the references for named set "a"; any read that maps to foo.fasta or bar.fasta will be considered a member of set a.
out_<name>=<output file>	Sets the output file name for reads mapping to set <name>.  out_a=stuff.sam would capture all the reads mapping to ref_a.
basename=<example%.sam>		This shorthand for mass-specifying all output files, where the % symbol is a wildcard for the set name.  For example, "ref_a=a.fa ref_b=b.fa basename=mapped_%.sam" would expand to "ref_a=a.fa ref_b=b.fa out_a=mapped_a.sam out_b=mapped_b.sam"
ref=<fasta files>		When run through the splitter, this is shorthand for writing a bunch of ref_<name> entries.  "ref=a.fa,b.fa" would expand to "ref_a=a.fa ref_b=b.fa".


Formats and Extensions
.gz,.gzip,.zip,.bz2	These file extensions are allowed on input and output files and will force reading/writing compressed data.
.fa,.fasta,.txt,.fq,.fastq	These file extensions are allowed on input and output files.  Having one is REQUIRED.  So, reads.fq and reads.fq.zip are valid, but reads.zip is NOT valid.  Note that outputting in fasta or fastq will not retain mapping locations.
.sam			This is only allowed on output files.
.bam			This is allowed on output files if samtools is installed.  Beware of memory usage; samtools will run in a subprocess, and it can consume over 1kb per scaffold of the reference genome.


Different versions:
BBMap			(bbmap.sh) Fastest version.  Finds single best mapping location.
BBMapPacBio		(mapPacBio.sh) Optimized for PacBio's error profile (more indels, fewer substitutions).  Finds single best mapping location.  PacBio reads should be in fasta format.
BBMapPacBioSkimmer	(bbmapskimmer.sh) Designed to find ALL mapping locations with alignment score above a certain threshold; also optimized for Pac Bio reads.
BBSplitter		(bbsplit.sh) Uses BBMap or BBMapPacBio to map to multiple references simultaneously, and output the reads to the file corresponding to the best-matching reference.  Designed to split metagenomes or contaminated datasets prior to assembly.
BBWrap			(bbwrap.sh) Maps multiple read files to the same reference, producing one sam file per input file.  The advantage is that the reference/index only needs to be read once.




Notes.
File types are autodetected by parsing the filename.  So you can name files, say, out.fq.gz or out.fastq.gz or reads1.fasta.bz2 or data.sam and it will work as long as the extensions are correct.

