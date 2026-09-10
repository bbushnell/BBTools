#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 11, 2026

Description:  Generates synthetic reads from a set of fasta assemblies.
Each assembly is assigned a random coverage level, with optional custom 
coverage for specific genomes.  Reads headers will contain the TaxID
of the originating genome, if the filename starts with 'tid_x_',
where x is a positive integer.

Default header style, where all numbers are 0-based:
f_(file number in argument list)
c_(contig number in the file)
s_(strand, 0 for + and 1 for -)
p_(start position on contig)
i_(insert size, useful for paired reads)
r_(reference length accounting for indels)
d_(1 if a PCR duplicate)
tid_(taxID, if present in the file or contig name)

Usage:  randomreadsmg.sh *.fa out=reads.fq.gz
or
randomreadsmg.sh ecoli.fa=40 mruber.fa=0.1 phix.fa=10 out=reads.fq.gz

File parameters:
in=<file,file>  Assembly input.  Can be a single file, a directory of files,
                or comma-delimited list.  Unrecognized arguments with no '='
                sign will also be treated as input files.
out=<file>      Synthetic read output destination.
out2=<file>     Read 2 output if twin files are desired for paired reads.

Processing parameters:
mindepth=1      Minimum assembly average depth.
maxdepth=256    Maximum assembly average depth.
depth=          Sets minimum and maximum to the same level.
reads=-1        If positive, set depth based on read length and genome size,
                to yield approximately this number of reads per file.
		Requires reading the input twice.
readspercontig=-1    If positive, ignore depth and make this many reads per contig.
mode=min4       Random depth distribution; can be min4, exp, root, or uniform.
cov_x=          Set a custom coverage level for the file named x.
                x can alternatively be the taxID if the filename starts
                with tid_x_; e.g. cov_foo.fa=5 for foo.fa, or cov_7=5
                for file tid_7_foo.fa
<file>=x        Alternate way to set custom depth; file will get depth x.
circular=f      Treat each contig as circular, and create spanning reads.
threads=        Set the max number of threads; default is logical core count.
                By default each input file uses 1 thread.  This flag will
                also force multithreaded processing when there is exactly 1
                input file, increasing speed for a complex simulation.
seed=-1         If non-negative, use the specified RNG seed.  Output content is
                then deterministic regardless of thread count: each input file's
                reads are a pure function of (seed, seed2, filename).

Multi-sample simulation parameters:
These allow generating multiple correlated samples from one genome set, for
testing differential-coverage binners.  Each sample is one invocation.
depthseed=-1    Seed for per-file base depth assignment; -1 uses 'seed'.
                Invocations sharing a depthseed (with different seeds) draw the
                SAME base depth per genome, so samples are correlated.
jitter=0.0      Multiply each file's depth by a random factor, symmetric in log
                space; jitter=0.1 gives roughly +-10% per-sample variation
                around the base depth.  With jitter=0 and a shared depthseed,
                samples are exact-duplicate in depth (fully correlated).
seed2=-1        Optional separate seed for the jitter stream; -1 uses 'seed'.
                Only needed to reproduce one specific jitter pattern
                independently of the generation seed.
zeroprob=0      Probability that each genome is ABSENT (depth 0) from this
                sample.  Keyed on depthseed, so samples sharing a depthseed
                share their presence/absence pattern while independent
                depthseeds draw independently - mimics sparse communities
                (e.g. NEON soil) that are mostly zeros in every library.
                Applies only to randomly-chosen depths, not custom or reads=.
Example - four correlated samples, ~10% depth wiggle, same community:
  for S in 1 2 3 4; do
    randomreadsmg.sh config=genomes.txt depthseed=777 seed=\$S jitter=0.1 \\
      out=sample\$S.fq.gz
  done

Artifact parameters
pcr=0.0         Add PCR duplicates at this rate (0-1).
randomkmer=f    Bias read start sites with random kmer priming.
kprime=6        Length for random kmer priming.
kpower=0.5      Raise linear primer distribution to this power (>0).
                Higher powers increase priming bias.
minkprob=0.1    Minimum primer kmer probability.

Platform parameters
illumina        Use Illumina length and error mode (default).
pacbio          Use PacBio HiFi length and error mode.
ont             Use ONT length and error mode.
paired=true     Generate paired reads in Illumina mode.
length=150      Read length; default is 150 for Illumina mode.
avginsert=300   Average insert size; only affects paired reads.

Long read parameters
minlen=1000     Minimum read length for PacBio/ONT modes.
meanlen=15000   Mean read length for PacBio/ONT modes.
maxlen=100000   Max read length for PacBio/ONT modes.
tailfactor=0.2  Controls heavy tail for ONT length distribution.
pbsigma=0.5     Log-normal standard deviation for PacBio length distribution.

Error parameters (all platforms)
adderrors=f     Set to true to add model-specific errors.
subrate=0.0     Add substitutions at this rate, independent of platform models.
insrate=0.0     Add length-1 insertions at this rate, independent of platform models.
delrate=0.0     Add length-1 deletions at this rate, independent of platform models.
indelrate=      Set insrate and delrate to half of this value.

Illumina-specific parameters
illuminanames=f Generate Illumina-format headers.
qavg=25         Average quality score, for generating Illumina errors.
qrange=0        Quality score range (+/- this much).
qflat=f         Use constant quality within a read, to increase compression
                when qrange>0.
addadapters=f   Add adapter sequence to paired reads with insert
                size shorter than read length.
adapter1=       Optionally specify a custom R1 adapter (as observed in R1).
adapter2=       Optionally specify a custom R2 adapter (as observed in R2).
illuminanames=f Make headers look like normal Illumina headers.
barcode=        Specify the barcode for Illumina headers.
machine=        Specify the machine for Illumina headers.

Long-read error parameters
Note: These may be overriden for any platform, including Illumina.
They are independent of, and applied in addition to, subrate/insrate/delrate.
srate=-1        Substitution rate; default 0.0025 ONT / 0.00015 PB.
irate=-1        Insertion rate; default 0.0055 ONT / 0.000055 PB.
drate=-1        Deletion rate; default 0.0045 ONT / 0.000045 PB.
hrate=-1        Homopolymer error boost; default 0.02 ONT / 0.000015 PB.
                The indel chance increases this much per homopolymer base.

Coverage variation parameters (used with 'sinewave' flag):
sinewave=f      Enable realistic coverage variation within contigs.
waves=4         Number of sine waves to combine; more waves create more 
                complex coverage patterns with irregular peaks and valleys.
waveamp=0.70    Controls the maximum variation in coverage due to the sine 
                waves.  Higher values (0-1) create more dramatic differences 
                between high and low coverage regions.
oribias=0.25    Strength of the origin of replication bias. Controls the max
                linear decrease in coverage from start to end of contigs.
minprob=0.10    Sets the minimum coverage probability as a fraction of target.
                Makes it improbable for regions have coverage that drops 
                below this level, preventing assembly gaps.
minperiod=2k    Minimum sine wave period, in bp.
maxperiod=80k   Maximum sine wave period, in bp.
variance=0.0    Vary coverage on a per-contig basis, within an assembly, by
                plus/minus this factor.  Unrelated to sinewave mode, which
		varies coverage WITHIN a contig.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

	parseJavaArgs "--xmx=1000m" "--xms=1000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP synth.RandomReadsMG $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" synth.RandomReadsMG "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
