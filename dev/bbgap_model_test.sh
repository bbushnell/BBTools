#!/usr/bin/env bash
# Full-reference API parity against frozen calls; run on Dori, not the laptop.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]]; then
	echo 'Usage: bbgap_model_test.sh models=PACKAGE in=QUERIES expected=CALLS count=1000 t=4 out=NEW_DIRECTORY'
	exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=256m --mode=fixed -eoom
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=4 -cp "$DIR/current" test.prot.ProteinFamilyAssignerTest "$@"
