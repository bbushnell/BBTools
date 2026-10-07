#!/bin/bash
# Inventories declared translation tables; never guesses a missing declaration.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: translationtablecensus.sh indir=paired_genomes out=census.tsv t=8 -Xmx1g'
  echo 'Scans a flat directory of .gff/.gff.gz files. Table 0 means missing metadata.'
  echo 'Reports CDS-row counts, missing/mixed codes, exception rows and exact FASTA pairs.'
  echo 'Does not validate biological assignments or select training examples.'
  exit 0
fi
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current" prok.TranslationTableCensus "$@"
