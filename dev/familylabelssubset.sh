#!/usr/bin/env bash
# Transfers labels to verified unchanged ordered model subsets; never infers new names.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
	cat <<'HELP'
Usage: familylabelssubset.sh parent=PARENT_MODELS models=SUBSET_MODELS
         labels=PARENT_LABELS out=NEW_RAW_TSV
Validates both model packages and matching protein sequences, preserves all
selected label fields, and rebinds local indexes and identity/consensus pins.
Output must be a new uncompressed TSV outside either model directory.
Run on the cluster for large references. Default heap: 8g; one worker.
HELP
	exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case ${arg,,} in
		--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
		-ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
		simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
		*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prok.FamilyLabelsSubset "${ARGS[@]}"
