#!/usr/bin/env bash
# Development experiment only; not a production model-reduction stage.
set -eo pipefail
if [[ $# -eq 0 || "$1" == --help || "$1" == -h ]]; then
	cat <<'HELP'
Development-only sequential consensus regeneration with member-support protection.
Uses public trnaconsensus.sh and the existing exhaustive Quantum support engine.

in=genes.fa ref=consensus.fa hbm=models.hbmc out=fresh_directory
expectedgenes=N expectedmodels=N cutoff=0.70 t=16
consensusconfig=single_group.config
graph=members.tsv       Optional saved support graph for the exact input library/corpus.
                       Without one, the existing engine performs the initial scan.
childheapmib=1024       Heap MiB for each serial public consensus invocation.

The consensus recipe must retain a single group (cluster=f, mincluster=1,
mingroup=1, prune=0) and keep its scientific settings explicit.
Model order and count stay fixed. Empty models remain unchanged. Each candidate
is accepted only if every former member still has support in the current library.
One final whole-corpus scan reports migrations and checks support retention.
Outputs include consensus.fa, models.hbmc, per-model decisions and assignments.
This experiment does not run a genomic caller or establish accuracy by itself.

mode=protect performs ONE restoration step from saved ordinary caller GFFs:
in=current.fa hbm=current.hbmc ref=original.fa originalhbm=original.hbmc
baseline=original.calls.gff calls=current.calls.gff truth=truth.tsv
family=SAM-I truthfamily=RF00162 expectedmodels=N out=fresh_directory
Lost references restore their baseline winner's original paired model. Other
pairs remain current. Output metrics describe the supplied calls, not the proposed
library. Run the ordinary caller on the proposal before deciding whether losses
are gone. Remaining losses with no newly revertible winner fail loudly.
Optional genes=corpus.fa graph=original.members.tsv expectedgenes=N cutoff=0.70
requests one final existing-engine support scan ONLY when the supplied calls have
zero baseline losses. support.summary.tsv reports named corpus losses/gains and
retention_pass separately; zero caller losses alone does not establish acceptance.
HELP
	exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" \
	prok.NcrnaSequentialRegenerationExperiment installation="$DIR" "$@"
