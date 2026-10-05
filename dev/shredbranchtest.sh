#!/bin/bash
# Independent branch oracle and exact real-genome partition/length audit.
set -o pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd) || exit 1
if [[ -d "$DIR/current" ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh" || exit 1
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case "${arg,,}" in
	--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx[0-9]*|xmx[0-9]*|-xmx=*|xmx=*|\
	-xms[0-9]*|xms[0-9]*|-xms=*|xms=*|\
	-ea|-da|ea|da|-eoom|eoom|-exitonoutofmemoryerror|exitonoutofmemoryerror|\
	simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("${arg,,}");;
	*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "${JVM_ARGS[@]}" || exit 1
setEnvironment || exit 1
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" synth.ShredBranchTest "${ARGS[@]}"
