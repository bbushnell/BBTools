#!/bin/bash
# Native SubnetGpuRoundTripTest validation with standard BBTools JVM setup.
set -o pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")" && pwd) || exit 1
if [[ ! -f "$DIR/javasetup.sh" ]]; then DIR=$(dirname "$DIR"); fi
if [ -f "$DIR/bbtools.jar" ]; then CP=$DIR/bbtools.jar; else CP=$DIR/current/; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs "--xmx=512m" "--xms=16m" "--mode=fixed" "$@"
setEnvironment
ARGS=()
CLASS=test.prot.SubnetGpuRoundTripTest
for arg in "$@"; do
	case "$arg" in
	selftest=t) CLASS=test.prot.SubnetGpuRoundTripFixtureTest;;
	--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-[xX][mM][xX]*|[xX][mM][xX]=*|-[xX][mM][sS]*|[xX][mM][sS]=*|\
	-da|-ea|da|ea|ExitOnOutOfMemoryError|exitonoutofmemoryerror|eoom|-ExitOnOutOfMemoryError|-exitonoutofmemoryerror|-eoom) continue;;
	*) ARGS+=("$arg");;
	esac
done
java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
