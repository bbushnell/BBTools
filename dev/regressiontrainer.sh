#!/bin/bash
#Written by Amber
#July 4, 2026

usage(){
echo "
Written by Amber
Last modified July 4, 2026

Description:  Trains a feed-forward neural network for CONTINUOUS outputs
              (regression) and saves it in BBNet format, readable by
              CellNetParser and everything else that consumes BBNets.
              Unlike train.sh, this is regression-native: pure MSE loss,
              Adam optimizer with cosine decay, no class balancing, no
              cutoff/FPR machinery.  Input standardization is learned from
              the data and folded into the first layer at export, so the
              saved network takes raw inputs.  After writing, the network
              is reloaded and verified against the in-memory model
              (round-trip check printed to stderr).

Usage:  regressiontrainer.sh in=<data.tsv> out=<net.bbnet> dims=16,32,1

Input format (same as train.sh / seqtovec.sh):
  First line '#dims <inputs> <outputs>', then tab-delimited floats,
  inputs first, target last.  One output supported; targets should be
  in [0,1] for final=sig, or roughly [-1.7,1.7] reachable for final=rslog.

Parameters:
in=<file>       Tab-delimited training vectors.
out=<file>      Output network (BBNet format).
dims=16,32,1    Network dimensions; hidden layers are tanh.
final=rslog     Final-layer activation: rslog (default) or sig.
                rslog is strongly recommended when targets touch the
                range boundary (sigmoid starves gradients there).
epochs=60       Training epochs (full passes).
batch=8192      Mini-batch size.
lr=0.003        Adam learning rate (cosine-decayed to 0).
wd=1e-4         Weight decay.
seed=1          RNG seed (deterministic).
vfraction=0.1   Validation split; best-validation network is exported.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")/.."
DIR="$(pwd)/"
popd > /dev/null

CP="$DIR""current/"

z="-Xmx4g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

regressiontrainer() {
	local CMD="java $EA $EOOM $z -cp $CP ml.RegressionTrainer $@"
	echo $CMD >&2
	eval $CMD
}

regressiontrainer "$@"
