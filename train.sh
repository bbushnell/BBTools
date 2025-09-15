#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 15, 2025

Description:  Trains or evaluates neural networks.

Usage:  train.sh in=<data> dims=<X,Y,Z> out=<trained network>

train.sh in=<data> netin=<network> evaluate

Input may be fasta or fastq, compressed or uncompressed.


I/O parameters:
in=<file>       Tab-delimited data vectors.  The first line should look like
                '#dims	5	1' with the number of inputs and outputs; the
                first X columns are inputs, and the last Y the desired result.
                Subsequent lines are tab-delimited floating point numbers.
                Can be created via seqtovec.sh.
validate=<file> Optional validation dataset used exclusively for evaluation.
net=<file>      Optional input network to train.
out=<file>      Final output network after the last epoch.
outb=<file>     Best discovered network according to evaluation metrics.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
evaluate=f      Don't do any training, just evaluate the network.
dims=           Set network dimensions.  E.g. dims=5,12,7,1
mindims,maxdims These allow random dimensions, but the number of inputs and
                outputs must agree.  e.g. mindims=5,6,3,1 maxdims=5,18,15,1
batches=400k    Number of batches to train.
alpha=0.08      Amount to adjust weights during backpropagation.  Larger 
                numbers train faster but may not converge.
balance=0.2     If the positive and negative samples are unequal, make copies
                of whichever has fewer until this ratio is met.  1.0 would
                make an equal number of positive and negative samples.
density=1.0     Retain at least this fraction of edges.
edges=-1        If positive, cap the maximum number of edges.
dense=t         Set dense=f (or sparse) to process as a sparse network.
                Dense mode is fastest for fully- or mostly-connected networks;
                sparse becomes faster below 0.25 density or so.

Advanced training parameters
seed=-1         A positive seed will yield deterministic output;
                negative will use a random seed.  For multiple networks,
                each gets a different seed but you only need to set it once.
nets=1          Train this many networks concurrently (per cycle).  Only the
                best network will be reported, so training more networks will
                yield give a better result.  Higher increases memory use, but
                also can improve CPU utilization on many-threaded CPUs.
cycles=1        Each cycle trains 'nets' networks in parallel.
setsize=60000   Iterate through subsets of at most this size while training;
                larger makes batches take longer.
fpb=0.08        Only train this fraction of the subset per batch, prioritizing
                samples with the most error; larger is slower.

Evaluation parameters
vfraction=0.1   If no validation file is given, split off this fraction of the
                input dataset to use exclusively for validation.
inclusive=f     Use the full training dataset for validation.  Note that
                'evaluate' mode automatically used the input for validation.
cutoffeval=     Set the evaluation cutoff directly; any output above this
                cutoff will be considered positive, and below will be
                considered negative, when evaluating a sample.  This does not 
                affect training other than the printed results and the best 
                network selection.  Overrides fpr, fnr, and crossover.
crossover=1     Set 'cutoffeval' dynamically using the intersection of the
                FPR and FNR curves.  If false positives are 3x as detrimental
                as false negatives, set this at 3.0; if false negatives are 2x
                as bad as false positives, set this at 0.5, etc.
fpr=            Set 'cutoffeval' dynamically using this false positive rate.
fnr=            Set 'cutoffeval' dynamically using this false negative rate.

Activation functions; fractions are relative and don't need to add to 1.
sig=0.6         Fraction of nodes using sigmoid function.
tanh=0.4        Fraction of nodes using tanh function.
rslog=0.02      Fraction of nodes using rotationally symmetric log.
msig=0.02       Fraction of nodes using mirrored sigmoid.
swish=0.0       Fraction of nodes using swish.
esig=0.0        Fraction of nodes using extended sigmoid.
emsig=0.0       Fraction of nodes using extended mirrored sigmoid.
bell=0.0        Fraction of nodes using a bell curve.
max=0.0         Fraction of nodes using a max function (TODO).
final=rslog     Type of function used in the final layer.

Exotic parameters
scan=0          Test this many seeds initially before picking one to train.
scanbatches=1k  Evaluate scanned seeds at this point to choose the best.
simd=f          Use SIMD instructions for greater speed; requires Java 18+.
cutoffbackprop=0.5   Optimize around this point for separating positive and
                     negative results.  Unrelated to cutoffeval.
pem=1.0         Positive error mult; when value>target, multiply the error 
                by this number to adjust the backpropagation penalty.
nem=1.0         Negative error mult; when value<target, multiply the error 
                by this number to adjust the backpropagation penalty.
fpem=10.5       False positive error mult; when target<cutoffbackprop
                value>(cutoffbackprop-spread), multiply error by this.
fnem=10.5       False negative error mult; when target>cutoffbackprop
                value<(cutoffbackprop+spread), multiply error by this.
spread=0.05     Allows applying fnem/fpem multipliers to values that
                are barely onsides, but too close to the cutoff.
epem=0.2        Excess positive error mult; error multiplier when 
                target>cutoff and value>target (overshot the target).
enem=0.2        Error multiplier when target<cutoff and value<target.
epm=0.2         Excess pivot mult; lower numbers give less priority to
                training samples that are excessively positive or negative.
cutoff=         Set both cutoffbackprop and cutoffeval.
ptriage=0.0001  Ignore this fraction of positive samples as untrainable.
ntriage=0.0005  Ignore this fraction of negative samples as untrainable.
anneal=0.003    Randomize weights by this much to avoid local minimae.
annealprob=.225 Probability of any given weight being annealed per batch.
ebs=1           (edgeblocksize) 8x gives best performance with AVX256 in
                sparse networks.  4x may be useful for raw sequence. 

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

z="-Xmx8g"
z2="-Xms8g"
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
	freeRam 8000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

train() {
	local CMD="java $EA $EOOM $z $z2 $SIMD -cp $CP ml.Trainer $@"
	echo $CMD >&2
	eval $CMD
}

train "$@"
