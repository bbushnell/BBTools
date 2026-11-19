#!/bin/bash
# Literal full-sequence, corrected-coordinate and saved-stratum integration.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_locus_overlap_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks collision-safe normalized exact DNA, both query/capture orientations,
inclusive containment, recovered/clan coordinate correction, orphan exclusion,
cross-assembly exact sequence, explicit unknowns, CLI parity and generated labels
through the actual saved-rank reporter. Also checks kept-only observer callbacks
and preservation of runtime/assertion failures after native input draining.
All genome paths are nonexistent. No caller, genome or alignment runs.
Compile the complete canonical current/ tree on Dori before execution.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlLocusOverlapTest "$@"
