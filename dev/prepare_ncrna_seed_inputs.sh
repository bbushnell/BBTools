#!/bin/bash
# Original genes absent from saved construction sequences; metadata only.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/prepare_ncrna_seed_inputs.sh in=prepared_population_prefix out=fresh_prefix expectedgenes=N expectedframes=N expectedmembers=N expectedmissing=N t=1 -Xmx1g -eoom
Reads the prepared .frames.tsv and .original_calls.tsv.gz ledgers.
Outputs .requests.tsv, .sources.tsv and .summary.tsv without opening genome FASTAs.
Only excluded truncated originals request new evaluation-only DNA; eligible exact
copies retain their saved same-group representative. Component merges are rejected.
Source size and availability are saved metadata, not fresh DNA validation.
No extraction, taxonomy, seed generation, caller or coverage is performed.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.PrepareNcrnaSeedInputs "$@"
