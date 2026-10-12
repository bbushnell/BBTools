#!/usr/bin/env bash
# Adds ProkCC native-harness bindings to a bare explicit-input MAG-QC final net.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case "${arg,,}" in
	--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
	-ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
	simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
	*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.MagQCFinalNetBinder "${ARGS[@]}"
