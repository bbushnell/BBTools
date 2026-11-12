#!/bin/bash
usage(){
cat <<'HELP'
Usage: dev/ncrnasweepfrontier.sh in=searches.tsv taxonomy=full.selection.tsv out=report

Reports saved native search policies without caller execution or payload reads.
Search manifest header (tab-separated): group, family, search, root.
root= values in that table name saved search directories containing contract.txt,
identity.tsv, caller.config, reference.tsv, familymap.tsv and roundN/ tables.
Groups must have identical declarations, caller/reference/map bytes and relevant
caller settings. Incompatible populations/runtimes must use separate groups.
Each completed native round is read once; incomplete rounds appear in .gaps.tsv.
Missing optional work/CDS counters are NA, never estimated or treated as zero.

taxonomy=       Complete original-input selection ledger, including unselected
                inputs, with input_id, fasta, domain, phylum, taxonomy_status.
                Every original ID must occur once. Exact ID and FASTA join is
                required for every reported input. No taxa are inferred.
out=prefix      .all.tsv retains every measured policy with search/round origin.
                .frontier.tsv retains unique nondominated points and shipped.
                .strata.tsv splits first gains/last losses for EVERY unique
                measured policy, including policies dominated on the total axes.
                .taxonomy.tsv retains corpus frequencies and category per input.
                .groups.tsv records input/settings provenance; .gaps.tsv names
                incomplete rounds. Existing outputs are protected by ow=f.

TP frontier: maximize TP, minimize aligned model bases, one strict improvement.
Breadth frontier additionally maximizes first-gain inputs and minimizes last-loss
inputs relative to shipped. Apparent FP/CDS are shown, never filtering axes.
Identical policy tuples deduplicate only after their full reported per-input and
summary counters agree; conflicting tuples remain separate and clearly marked.
Equality of these counters does not prove equality of the underlying call loci.
The lower quartile uses nearest rank over distinct informative domain/phylum
labels, including all ties. UNKNOWN and UNCLASSIFIED have separate categories.
This measures representation in this corpus, not global sequencing rarity.

Logical alignments/model bases are work proxies. Model bases omit query length;
cached CPU does not measure whole-caller time. Whole-caller seconds stay UNMEASURED.
No numerical operating point is selected or adopted. One report thread; retained
primitive vectors scale with inputs times unique policies within one group.
HELP
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaSweepFrontier "$@"
