#!/usr/bin/env bash
# Select original whole consensus/HBM pairs by retained FASTA and added names.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/exportncrnamodelsubset.sh in=original.fa hbm=original.hbmc keep=retained.fa names=additions.txt out=fresh_directory
The retained FASTA must preserve original names, complete headers and bases.
The additions file contains one unique native model ID per line, disjoint from
retained names; it may be empty. Unknown or duplicate names fail before export.
Every original FASTA/HBM pair is validated, including unselected pairs. Existing
paired-record selection preserves original order, headers, upper-case sequences
and all positional counts. Written pairs are decoded and checked in full.
FASTA wrapping or HBM encoding can normalize; whole-file byte identity is not
claimed. This exports a paired library, not an installable family or a verified
caller result. No alignment, fitting, scoring or gene calling is performed.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ExportNcrnaModelSubset "$@"
