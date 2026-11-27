#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/bindordinalncrnacalls.sh in=saved.gff bindings=panel.sources.tsv callerconfig=actual.config expectedinputs=1000 out=fresh_directory
Splits saved GFF input_ordinal0 blocks into per-input query views without DNA or
reference reads. Binds ordinals through input_ordinal0/assembly_id/input_id/fasta/
canonical_new_gff metadata, retaining the exact global input ID and FASTA path.
The original flat caller configuration must list those inputs in the same order
and point out= to the saved GFF. Boundaries must be contiguous from zero.
Every feature row and its order is preserved; query headers are normalized and
source_lines.tsv maps each new feature line back to its original GFF line.
Set compressed=t to write queries/*.gff.gz and source_lines.tsv.gz through
native Java streams, without a compressor subprocess per input. The default
compressed=f preserves plain query files and the original source-line filename.
Small manifests remain plain and point explicitly to the selected query paths.
Compressed mode rejects forcebgzip=t and forcepigz=t before creating outputs.
Zero-feature inputs retain a file and manifest row. calls.tsv and references.tsv
feed GradeNcrnaReference; this binder assigns no families and performs no match.
The inputs and run directory identify the binding; no provenance-only hash is
required. Version-2 complete.tsv records source path and counts and is written
last; partial directories are failures. The removed inputsha80 option is rejected.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.BindOrdinalNcrnaCalls "$@"
