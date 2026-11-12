#!/bin/bash
# Saved-call projection, preserving explicit input ownership and physical lines.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 8, 2026

Usage: dev/bindncrnacalls.sh in=combined.gff.gz genomes=panel.tsv \
       panelcontigs=panel_contigs.tsv observed=contigs.tsv.gz census=genomes.tsv \
       memberships=memberships.tsv callerconfig=original.config baseline=name \
       inputsha80=20hex out=new_directory

Projects EVERY generic RNA row from a saved combined CallGenes GFF into per-input
query GFFs, including duplicates and unknown/excluded models. No new calls, genome
reads, model attribution, scoring or biological reference decisions. Other feature
types stay in the immutable original; their counts remain in the output inventory.

Required input schemas (tab-separated headers, shown with spaces):
genomes= index genome assembly taxid domain phylum gff fasta loci genome_bases contigs
         Original panel, including inputs with no RNA calls.
panelcontigs= genome assembly contig length
         The cobalamin spelling seqid in place of contig is also accepted.
         Literal contig names must be unique across the entire panel.
census=  input_id fasta file_bytes nucleotide_bases fasta_records
memberships= source member_index input_id supplied_fasta source_genome source_assembly taxid declared_compartment reported_gff reported_sha80
         Existing evidence catalog. Only exact canonical or recorded membership
         FASTA paths bind; metadata agreement is not independent sequence proof.
observed= input_id fasta literal_target observed_length
         Saved observations, evidence TSV escaping (%25/%09/%0D/%0A). Every panel
         contig must match its bound catalog input and observed length.
callerconfig= Flat original key=value file with one in= or in1= per panel input,
         in exact panel order. Nested config= and in2= inputs are rejected.
in=      Original nine-column CallGenes GFF; source must be BBTools. Every feature
         must have a known owner, valid positive bounds, strand +/-, and type CDS,
         tRNA, rRNA or RNA. Comments/blank lines are counted separately.
baseline= Explicit label identifying this saved caller configuration.
inputsha80= Declared 20-character lowercase fingerprint of the original GFF.
         Recorded only: the adapter does NOT compute or independently verify it.

Outputs (fresh directory only, no overwrite):
input_N.gff.gz      Generic RNA in original per-input order, unchanged row bytes
                   with LF newlines; one header, including for empty inputs.
calls.tsv          Literal UTF-8 input_id/fasta/gff query manifest, published LAST.
bindings.tsv       Original and canonical FASTA paths, panel source lines.
line_map.tsv.gz    Original GFF physical line, input_id, query physical line.
original_counts.tsv Per-input original CDS/tRNA/rRNA/RNA counts and total features.
contigs.tsv.gz     Exactly bound saved contig observations for this panel.
sources.tsv        Original metadata paths, baseline and declared GFF fingerprint.
summary.tsv        Global line/type counts, limits, peak retained array capacity.
                   Completion requires exit zero AND calls.tsv.

Limits (K/M/G suffixes accepted):
maxbuffer=134217728 Aggregate retained RNA byte-array and long-array CAPACITIES,
                   including spare capacity; checked before growth. Maps, native
                   I/O buffers and temporary allocation are separately heap-bound.
maxlines=2000000   Maximum retained RNA rows.
maxoutput=268435456 Maximum logical uncompressed bytes across ALL outputs.
                   This is not filesystem allocated space; gzip/filesystem overhead
                   is separate. Measure physical footprint for storage planning.
maxfiles=10000     Maximum number of output files: panel size plus seven. The
                   directory uses one additional inode. Failure can leave partial
                   files but never a completed query manifest.
Defaults: t=1, -Xmx1g, ow=f. For many small outputs, pigz=f bgzip=f avoids launching
one external compressor per input. Unicode paths require a supporting native locale.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.BindNcrnaCalls "$@"
