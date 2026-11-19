#!/usr/bin/env bash
# Convert canonical protein substitution matrices between native TSV and NCBI text.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: convertproteinmatrix.sh in=MATRIX out=OUTPUT outformat=ncbi|native

Input format is detected from content, independently of its extension.
Native version1 and self-described NCBI files retain their name/scale/gap.
Plain NCBI input requires name=IDENTIFIER scale=N gap=N; conflicting metadata fails.
NCBI output additionally requires background=PROBABILITIES.tsv band=NAME.
The background is a versioned normalized monomer probability table, using
(query_marginal+reference_marginal)/(2*total). No frequencies are guessed from scores.
POSSUM55 uses its accepted 25to55 band, scale=4 and gap=16. Its name expands to
Prokaryotic Ortholog Sequence SUbstitution Matrix; scores are round(8*log2(O/E))
in eighth-bit units. Native scale4 describes the multiplier relative to BLOSUM units.
The converter verifies the supplied25to55 table reproduces all400 scores before
adding that derivation description. This verifies the pairing, not its biological origin.

All400 canonical scores are retained exactly. Exported B=D/N and Z=E/Q are
background-weighted; X is the background expected score. Ambiguous pairs use
both distributions and one rounding, half away from zero. Star uses the canonical
matrix minimum against residues and +1 against itself. Recommended suffix: .mat.
Import validates B/Z/X/star/J cells but discards them from runtime binding.
BBTools still folds B/Z/J to X, uses scaled BLOSUM62 X scores, rejects internal
stops and strips one conventional terminal FASTA stop. Export does not change that.
Statistical E-value parameters are not provided. This is text interchange,
not a claim that BLAST supports this custom matrix or its statistical calibration.

Defaults: outformat=ncbi, overwrite=f, one conversion thread, -Xmx256m.
Use ow=t to replace an output; input/output aliases and append mode are rejected.
Native output needs no background. Both formats can use ordinary compression suffixes.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case ${arg,,} in
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount=1 -cp "$CP" prot.ProteinSubstitutionMatrixIO "${ARGS[@]}"
