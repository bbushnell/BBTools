#!/usr/bin/env bash
# Experimental position-specific profile scoring; production defaults are untouched.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Experimental HBM scoring; this tool never activates a ProkCC model or changes its defaults.
Usage: hbmpositionassay.sh mode=run in=proteins.faa pairs=frozen_top50.tsv resources=DIR out=PREFIX
  path=blosum|profile|both selects ordinary BLOSUM, profile-aware alignment, or both traces.
  Selected profile recipe: kind=logodds beta=0.01 clip=t clipmin=-4 path=profile
  The experimental defaults remain beta=0.1 clip=f path=both; pass the selected recipe explicitly.
  refs=FASTA hbm=MODEL provenance=TSV override the experimental input artifacts.
  background=TSV supplies the original frozen 20-residue background for rebuild comparisons.
  pairs= uses the saved ten-column shortlist with exactly 50 unique family IDs per query.
  t=1 sets workers; validate=t checks each profile traceback against its optimum.
  Outputs: PREFIX.tsv and PREFIX.time.tsv. Existing files are refused.
Other modes: mode=test; mode=grade root=REFERENCE_DIR in=RESULT.tsv out=PREFIX;
  mode=gradeaccepted root=REFERENCE_DIR in=RESULT.tsv out=PREFIX n=10000
Scores are experimental profile units, not calibrated production thresholds or BLOSUM E-values.
Production assignment and calibrated resource integration are deferred to ProkCC 1.5.
Legacy positional run, grade and test commands remain supported.
HELP
 exit 0
fi
kind=${1#mode=}; shift
case $kind in run|grade|gradeaccepted) main=prot.HbmPositionAssay;; test) main=prot.HbmPositionTest;; *) echo 'Expected mode=run|grade|gradeaccepted|test' >&2; exit 2;; esac
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
explicit_mode=false
for arg in "$@"; do
	case ${arg,,} in
	 --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
	 -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
	 simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
	 mode=*) explicit_mode=true; ARGS+=("$arg");;
	 *) ARGS+=("$arg");;
	esac
done
if [[ $explicit_mode == false && ( $kind == grade || $kind == gradeaccepted ) ]]; then ARGS=("mode=$kind" "${ARGS[@]}"); fi
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$main" "${ARGS[@]}"
