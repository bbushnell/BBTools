#!/bin/bash
# Metadata-only first-fit batches for unambiguous saved-call binding.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden

Usage: dev/partitionncrnabenchmarkinputs.sh in=panel.1000.selection.tsv \
       census=genomes.tsv memberships=memberships.tsv \
       contigs=evidence.contigs.tsv.gz out=new_directory

Keeps each selected original input intact, assigning it to the first batch in
which none of its literal contig names occurs. Census order defines deterministic
placement. Reads saved metadata only; never opens or renames a FASTA, GFF or model.
Use the SAME recorded partition for every caller build being compared.

in=           Eight-column selection table from dev/selectncrnabenchmarkinputs.sh:
              input_id fasta nucleotide_bases fasta_records domain phylum
              base_log2 record_log2
census=       Native five-column census: input_id fasta file_bytes
              nucleotide_bases fasta_records
memberships=  Native ten-column source-membership table from that same census.
contigs=      Complete saved observation table: input_id fasta literal_target
              observed_length. FASTA/name fields use native percent/tab/CR/LF
              escaping. Worker-buffered rows may interleave different inputs.
              Every catalog input must be present with exact path/base/record
              dimensions. Decoded names cannot contain tab, CR or LF.
out=          Fresh directory under an existing parent. No overwrite/append.

Outputs:
assignments.tsv       Every selected original ID/path, batch and local ordinal.
batches.tsv           Each batch's input/base/record totals and binder file budget.
batch_N.inputs.config Flat in= declarations in original census order.
batch_N.genomes.tsv   Eleven-column BindNcrnaCalls genome roster.
batch_N.contigs.tsv   Four-column literal-name/length roster for that binder.
batch_N.binder.config Explicit maxfiles=batch input count+7 for binder metadata.
summary.tsv           Input/base/record totals, distinct/shared names, shared-name
                      input/occurrence counts, multiplicity lower bound and
                      greedy batches.
                      Written only after all batch writers close successfully.

Genome labels are original catalog input IDs; assembly and unknown historical/
biological fields are NA. These are generated binding labels, not source claims.
Within selected inputs, duplicate names are rejected: no whole-input partition
can fix them. Unselected names are not interned or checked for duplication;
complete-catalog path, base and record dimensions are still checked.
Names are compared after decoding, with UTF-8 bytes preserved. The partition is
valid but not necessarily minimal. Maximum input multiplicity of any one name
is a lower bound; it need not equal the achieved first-fit batch count.

For each untimed grading batch, combine its inputs with the exact flat caller
settings, run stock callgenes.sh once, then use dev/bindncrnacalls.sh with the saved
GFF, generated genome/contig rosters, original observed/census/membership tables
and that exact flat caller config. Include the generated binder config so the
default 10,000-file budget does not exclude a large batch. The caller build and
original-GFF identity belong to the owning recipe. This tool performs no calling,
binding, grading, biological admission or input identity authentication.

Zero exit and NCRNA_BENCHMARK_PARTITION_COMPLETE are required. Partial outputs
may remain after an I/O error; an existing summary alone is not a success check.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.PartitionNcrnaBenchmarkInputs "$@"
