#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 6, 2026

Description:  Synthesizes MAG-QC (CheckM3) neural-net TRAINING vectors.  For each
synthetic bin it samples one target organism's contigs to a drawn completeness,
adds contaminant contigs to a drawn contamination, and emits the bin's
protein-family count vector plus its ACHIEVED (completeness, contamination)
regression targets.  Organisms are split into train/val BEFORE sampling so a
held-out organism never appears in the training file.  Reads a per-contig
precompute cache (built by the mmseqs pipeline), NOT the native aligner.

This is a DEVELOPMENT tool (training-data generation); it is not part of the
end-user MAG-QC path.

Usage:  magqcvectormaker.sh cache=<cache.tsv> sizemap=<sizemap.tsv> \\
          taxpgm=<taxid_pgm.tsv> out=<train.tsv> [outval=<val.tsv>] [options]

Required parameters:
cache=<file>    Per-contig precompute cache TSV.
sizemap=<file>  tid<tab>genome_bp map.
taxpgm=<file>   tid<tab>phylum<tab>pgm map (defines the phylum one-hot).
out=<file>      Output training-vector TSV (has a '#dims' header line).

Optional parameters (and their defaults):
outval=<file>   Output held-out-organism validation TSV.
familylist=<f>  Family list (sets the number of family columns).
features=<file> Reduced-feature rank file (one kept family rank per line); emits
                only those family columns.
tree=<file>     TaxTree (enables same-family contaminant bias).
n=400000        Training bins to generate.
valn=40000      Validation bins to generate.
valfrac=0.10    Fraction of organisms held out for validation.
seed=1          RNG seed.
enc=ratio       Family-count encoding: ratio N/(1+N); raw min(N,32)/32;
                log log2(1+N)/log2(65); two 2-cols (presence+excess); norm
                0-if-absent-else N/avgCopyWhenPresent (per-family baseline).
minlen=0        Drop contigs shorter than this.
cleanspike=0.15 Probability a bin is drawn with zero contamination.
multicontamprob=0.15  Probability of 2-3 contaminants instead of 1.
samefamprob=0.70      Probability a contaminant is from the target's family.

Subnet emission (composite architecture; emits an ALIGNED per-bin subnet file
alongside the global vectors, same bins, global output byte-identical):
subnet=<name>   ncrna (observed ncRNA counts) or famset (observed counts of an
                arbitrary family-rank subset).  Target = the organism's NATIVE
                subset complement (the denominator; completeness is derived
                downstream as observed/predicted, see subnetratioscore.sh).
subnetout=<f>   Subnet training-row output TSV.
subnetvalout=<f>  Subnet validation-row output TSV.
subsetfile=<f>  famset only: the subset definition, one family rank per line
                (e.g. a per-phylum marker set or a co-occurrence module).
snbinscaled=f   Context: reparameterize bin size as log2(1+bp/2048)/16.
sndomain=f      Context: add an 8-way domain one-hot.
snhhcaga=f      Context: add per-organism HH/CAGA dimer features (needs kmerfile=).
sngenelen=f     Context: add mean gene length + gene-length stddev.
sncodingaffine=f  Context: affine-rescale coding density (no-op without wd).
kmerfile=<f>    tid<TAB>HH<TAB>CAGA table from genomedimerfeatures.sh.

Java Parameters:
-Xmx            Set Java heap (overrides autodetection); the cache is held in RAM.
-eoom           Exit if an out-of-memory exception occurs.
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=8g" "--xms=1g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.MagQCVectorMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
