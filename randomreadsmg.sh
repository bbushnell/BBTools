#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 15, 2025

Description:  Generates synthetic reads from a set of fasta assemblies.
Each assembly is assigned a random coverage level, with optional custom 
coverage for specific genomes.  Reads headers will contain the TaxID
of the originating genome, if the filename starts with 'tid_x_',
where x is a positive integer.

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
reads=-1        If positive, ignore depth and make this many reads per contig.
variance=0.5    Coverage within an assembly will vary by up to this much;
                one region can be up to this fraction deeper than another.
mode=min4       Random depth distribution; can be min4, exp, root, or linear.
cov_x=          Set a custom coverage level for the file named x.
                x can alternatively be the taxID if the filename starts
                with tid_x_; e.g. cov_foo.fa=5 for foo.fa, or cov_7=5
                for file tid_7_foo.fa
<file>=x        Alternate way to set custom depth; file will get depth x.
threads=        Set the max number of threads; default is logical core count.
                By default each input file uses 1 thread.  This flag will
                also force multithreaded processing when there is exactly 1
                input file, increasing speed for a complex simulation.
seed=-1         If positive, use the specified RNG seed.  This will cause
                deterministic output if threads=1.

Artifact parameters
pcr=0.0         Add PCR duplicates at this rate (0-1).
randomkmer=f    Bias read start sites with random kmer priming.
kprime=6        Length for random kmer priming.
kpower=0.5      Raise linear primer distribution to this power (>0).
                Higher powers increase priming bias.
minkprob=0.1    Minimum primer kmer probability.

Platform parameters
illumina        Use Illumina length and error mode (default).
pacbio          Use PacBio length and error mode.
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
indelrate=0.0   Add length-1 indels at this rate, independent of platform models.

Illumina-specific parameters
qavg=25         Average quality score, for generating Illumina errors.
qrange=0        Quality score range (+/- this much).
addadapters     Add adapter sequence to paired reads with insert
                size shorter than read length.
adapter1=       Optionally specify a custom R1 adapter (as observed in R1).
adapter2=       Optionally specify a custom R2 adapter (as observed in R2).

Long-read error parameters
Note: These may be overriden for any platform, including Illumina.
srate=-1        Substitution rate; default 0.0025 ONT / 0.00015 PB.
irate=-1        Insertion rate; default 0.0055 ONT / 0.000055 PB.
drate=-1        Deletion rate; default 0.0045 ONT / 0.000045 PB.
hrate=-1        Homopolymer error boost; default 0.02 ONT / 0.000015 PB.
                The indel chance increases this much per homopolymer base.

Coverage variation parameters (only used with 'sinewave' flag):
sinewave        Enable realistic coverage variation within contigs.
numwaves=4      Number of sine waves to combine; more waves create more 
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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    # Parse Java arguments with tool-specific defaults
    # Use auto mode with 84% of available RAM, minimum 1000MB
    parseJavaArgs "--mem=1000m" "--percent=84" "--mode=auto" "$@"
    
    # Set environment paths
    setEnvironment
    
    # Set the Java memory parameters
    z="-Xmx${RAM}m"
    z2="-Xms${RAM}m"
}
calcXmx "$@"

generate() {
	local CMD="java $EA $EOOM $XMX $XMS -cp $CP synth.RandomReadsMG $@"
	echo $CMD >&2
	eval $CMD
}

generate "$@"
