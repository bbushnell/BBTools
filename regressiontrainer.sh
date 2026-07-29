#!/bin/bash

usage(){
echo "
Written by Brian Bushnell, Amber, Neptune
Last modified July 29, 2026

Description:  Trains a feed-forward network for CONTINUOUS outputs and writes it in
BBNet format, byte-compatible with everything that consumes BBNets.  Differs from
train.sh (ml.Trainer): pure MSE regression with no balancing, cutoffs or FPR/FNR
machinery; Adam with cosine learning-rate decay; mini-batches; and input
standardization learned from the data and folded into the first layer at export, so
the saved net consumes RAW inputs.  Do not standardize upstream.

Input is streamed, so memory scales with the sample count rather than the file size.
After writing, the net is reloaded and checked against the in-memory model; the
'round-trip check' line must read (OK) or the output is not trustworthy.

Usage:  regressiontrainer.sh in=<data.tsv> out=<net.bbnet> dims=<16,32,1>

Example:
  regressiontrainer.sh in=training.tsv out=net.bbnet dims=48,64,32,1 epochs=100 final=sigmoid

Input format:  same as train.sh — a '#dims <in> <out>' header line, then
tab-delimited floats: <in> input columns followed by <out> target columns.

Required parameters:
in=<file>       Training vectors, tab-delimited with a '#dims' header.
out=<file>      Output network, in BBNet format.
dims=48,64,32,1 Layer widths: inputs, hidden..., outputs.  The last entry sets the
                output width; more than one output is supported and the loss is
                summed squared error across them.  Optional if netin= is given.

Training parameters:
epochs=60       Passes over the training data.
batch=8192      Mini-batch size.
lr=0.003        (alpha) Learning rate.  Set this explicitly; alpha is an alias and a
                stale alpha= in a script will silently override.
wd=1e-4         Weight decay.
seed=1          RNG seed; identical seeds reproduce identical nets.
vfraction=0.1   Fraction of data held out for validation.

Network shape:
final=rslog     Output activation: linear, rslog, or sigmoid.  DEFAULT IS RSLOG,
                WHICH IS UNBOUNDED.  Set it explicitly whenever the output range
                matters (for example when computing calibration error), or you may
                get values well outside [0,1] that look plausible.  Sigmoid requires
                targets in [0,1] and saturates at exactly 0 and 1.
hidden=tanh     Hidden activations.  One name is uniform; a comma-separated list
                assigns one per cell at random from the set, seeded and reproducible
                (e.g. hidden=tanh,swish,sig).  Names: tanh sig rslog msig swish esig
                emsig bell linear.  Default tanh is also the fastest path.
netin=<file>    Continue training from an existing .bbnet instead of random
                initialization.  dims= may then be omitted; the net supplies it.
                Training runs in raw input space because the saved net already has
                its standardization folded in.  The extraction is verified against
                the source net before training starts and refuses to run if it
                disagrees.

Performance:
simd=f          Train in float using vectorized kernels.  Measured 1.46x on a narrow
                net and 3.39x on a wide one, with no quality difference; benefit
                scales with layer width.  Results are close to but not bit-identical
                to the scalar path, so it is off by default.
pad8=f          Round hidden layers up to a multiple of 8.  Not a speedup: wall time
                is unchanged while the net does about 11% more arithmetic, so it buys
                capacity rather than time.

Options that were measured NEUTRAL or HARMFUL on synthetic data.  All default off.
Treat them as unproven rather than recommended, and measure on your own data:
density=1       Fraction of hidden-layer edges to keep; the output layer stays dense.
                Produces a smaller model but NOT a faster run, since pruned weights
                remain zeros in dense arrays.
density1=0      Density override for the first hidden layer; 0 uses density.
edgeblock=1     Block size for structured sparsity.  Use 8 to let sparse kernels
                vectorize.
sparse=f        Gather only live edges instead of multiplying through pruned zeros.
                Measured 26% SLOWER at density=0.75 and only 9% faster at
                density=0.25 on a wide net; the gather costs about what it saves.
sort=f          Prioritize high-error and stale samples instead of visiting all
                uniformly.  Measured a wash at equal compute.
setsize=0       Training samples visited per epoch when sort=t; 0 means all.
normalize=f     Z-score weight standardization during training.  Measured 88x WORSE
                on synthetic data at the default strength.
normfactor=0.125  Blend strength for normalize.
normshrink=0.999  Per-application decay of that strength.

Robustness:
checkexp=f      Scan input for scientific notation, which the fast float parser
                misreads silently.  Costs about 8% of load time, so it is off by
                default; BBTools tools do not emit exponent notation.  Turn it on for
                vector files produced elsewhere, or check yourself with
                'grep -c \"[eE]\" file.tsv' (must be 0).
fjpd=f          (forcejavaparsedouble) Force exact float parsing everywhere.
overwrite=t     (ow) Overwrite an existing output file.

Java Parameters:
-Xmx            Set max memory.  Vectors are held in RAM: budget roughly
                rows * inputs * 4 bytes, plus about 100 MB overhead.  Five million
                rows of 48 inputs needs about 1 GB, so -Xmx2g is comfortable.  Do not
                size from observed RSS under a default heap; the JVM takes far more
                than it needs and you will over-request several-fold.
-eoom           Exit on out-of-memory.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--percent=60" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.RegressionTrainer $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
