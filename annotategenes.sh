#!/usr/bin/env bash
# Existing gene calls gain HBM-family annotations; coordinates remain unchanged.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
	cat <<'HELP'
Usage: annotategenes.sh genome.fna calls.gff out=annotated.gff
   or: annotategenes.sh in=genome.fna gff=calls.gff out=annotated.gff

Annotates contiguous CDS calls using the existing expanded-family HBM package.
Reads nucleotide FASTA and GFF with BBTools readers; threads perform alignments.
Install the BBGAP reference tar by extracting it into BBTools/resources/.
It creates resources/bbgap/ with models, labels, evidence and usage instructions.
FASTA and GFF positional inputs are recognized by extension, in either order;
compressed extensions are supported. Each input may be specified only once.

models=DIR        Default: ?bbgap/models (searched in BBTools/resources/).
                  Keep the complete nine-file package and artifacts.tsv together.
labels=FILE       Default: labels/family_labels.tsv.gz beside the selected models
                  directory; uncompressed .tsv is also accepted. Labels must match
                  the model roster. labels=null disables names. Explicit models
                  without a sibling label file produce family assignments only.
                  Explicit models/labels paths also support ? resource lookup.
evidence=FILE     Optional per-CDS TSV with scores, coverage, table and status.
transl_table=11   Fallback table: 4, 11 or 25. Row transl_table overrides a file
                  ##transl_table directive, which overrides this fallback.
t=N              Alignment workers; default available processors.
ow=t             Overwrite outputs. Append and paired inputs are not supported.

Preserves GFF record order, comments, RNA features, CDS coordinates and original
attributes, including partial and anticodon. Appends bbgap_* attributes without
overwriting product. Empty attribute '.' becomes the added fields; line endings
are normalized to LF. Inputs already containing bbgap_* fields are rejected.
Output paths must be outside the immutable model package.

Statuses are assigned, no_hit (searched) or unsupported (not searched). Compound
CDSs, unknown phase/strand, unsupported tables, recoding exceptions, ambiguous
codons, internal stops and inconsistent CDS boundaries are unsupported. Invalid
GFF structure, missing contigs and out-of-bounds intervals fail loudly.

Binary partial flags refer to genomic left/right edges. Legacy partial=true
leaves both boundaries unknown. Complete starts use table-specific initiation M;
complete ends must contain a terminal stop, which is omitted from the peptide.
Three-prime partials retain only whole codons; a stop there is unsupported.

Search uses all top 50 F4 candidates, positional HBM scoring, identity >=40.612846%
and mutual paired-core coverage >=0.8, with deterministic representative-ID ties.
Partial genes use the same coverage gate. Scores are not E-values or calibrated
function probabilities. Evidence reference endpoints are zero-based inclusive.
Public database searches are used to build labels; this tool needs no web service.
HELP
	exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case ${arg,,} in
		--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
		-ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
		simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
		*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$DIR/current" prok.AnnotateGenes "${ARGS[@]}"
