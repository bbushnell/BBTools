#!/bin/bash
# Model-bound profile round-trip, score/order parity and rejection fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_profile_resource_test.sh out=fresh_fixture_directory t=1 -Xmx512m -eoom
   or: dev/ncrna_ddl_profile_resource_test.sh verify=tiny_builder_output.tsv t=1 -Xmx512m -eoom
Checks complete profile key/vote round-trip, actual raw scores and hybrid order,
independent worker scratch, UTF-8 model IDs, consensus content/order and malformed
metadata/count/key/row rejection. Empty and maximum16bit keys are literal tests.
No genome, verifier, caller or performance measurement.
The verify mode checks the actual builder CLI output from the prior accepted
two-member saved-window fixture against its literal model/member definitions.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlProfileResourceTest "$@"
