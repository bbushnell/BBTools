#!/usr/bin/env bash
# Exact union-versus-partition integration test, with one alignment worker.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
	cat <<'HELP'
Usage: proteinfamilycompetitiontest.sh models=UNION_MODELS packs=PART_A_MODELS,PART_B_MODELS
         in=proteins.faa count=N out=NEW_PASS_FILE
Partitions must preserve the union's complete family order and scientific inputs.
Checks legacy single-library, combined single-pack and combined multi-pack results.
Run large model tests on the cluster. Default heap: 8g; one worker thread.
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.ProteinFamilyCompetitionTest "${ARGS[@]}"
