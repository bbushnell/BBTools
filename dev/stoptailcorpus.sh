#!/usr/bin/env bash
# One-genome tail extraction, or finite literal regression fixtures. No caller jobs.
set -e
usage(){
 printf '%s\n' 'Usage: bash dev/stoptailcorpus.sh in=genome.fa gff=genes.gff out=tails.fa outmeta=tails.tsv outstats=stats.tsv genome=ID source=STRATUM taxid=N phylumtaxid=N taxnamespace=NAME split=build|calibration|final'
 printf '%s\n' 'Single-genome options: table=11 tableassumed=t/f minlen=100 cap=200 outaudit=cds_audit.tsv overwrite=f'
 printf '%s\n' 'Batch usage: bash dev/stoptailcorpus.sh mode=batch manifest=/path/genomes.tsv out=/path/NEW_DIR'
 printf '%s\n' 'Batch requires an absolute regular-file manifest and a fresh output directory. Each row supplies one genome; FASTA/GFF paths must be absolute regular files.'
 printf '%s\n' 'Batch TSV header/columns, in order: genome_id, source, fasta, gff, taxid, phylum_taxid, taxonomy_namespace, split, table, table_assumed.'
 printf '%s\n' 'Batch defaults: minlen=40 cap=200 maxoriginal=99 rejectpseudo=f reject5primepartial=f. Length and rejection options may be overridden explicitly.'
 printf '%s\n' 'Population filters: rejectpseudo=f reject5primepartial=f (also accepted in batch mode). Set rejectpseudo=t for both regenerated views; reject5primepartial=t for SHORT only.'
 printf '%s\n' 'Counters are exclusive outcomes: structural exclusions, pseudo policy, existing extraction checks, then five-prime policy. Three-prime partial CDS remain excluded in both views.'
 printf '%s\n' 'maxoriginal=N limits original annotation length after longest-eligible same-stop dedup; single-genome default unlimited, batch default99. It does not truncate long proteins.'
 printf '%s\n' 'Table is a fallback; feature/file declarations take precedence. MAG caller assumption flags remain assumed.'
 printf '%s\n' 'NCBI partial ranges and genomic edge bits are strand-aware. Only retained tails are validated; no forced initiator M.'
 printf '%s\n' 'Development-module-only tests: bash dev/stoptailcorpus.sh test NEW_FIXTURE_DIRECTORY [-Xmx512m]'
 printf '%s\n' 'Requires a separately supplied test module on the classpath; test.prok.StopTailCorpusTest is not included in the shipped selection.'
}
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then usage; exit 0; fi
script=$(readlink -f "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$script")/.." && pwd)
suite=prok.StopTailCorpus
fixture=
if [[ $1 == test ]]; then
 shift; fixture=${1:?new fixture directory}; shift
 suite=test.prok.StopTailCorpusTest
fi
module= dependency=
args=()
for arg in "$@"; do
 case "$arg" in
  mode=batch) suite=prok.StopTailCorpusBatch;;
  library=*) DIR=${arg#*=};;
  modulejar=*) module=${arg#*=};;
  dependencyjar=*) dependency=${arg#*=};;
  *) args+=("$arg");;
 esac
done
set -- "${args[@]}"
CP=$DIR/current/
[[ ! -f $DIR/bbtools.jar ]] || CP=$DIR/bbtools.jar
if [[ -n $dependency ]]; then [[ -f $dependency ]] || exit 2; CP=$dependency:$CP; fi
if [[ -n $module ]]; then [[ -f $module ]] || exit 2; CP=$module:$CP; fi
. "$DIR/javasetup.sh"
. "$DIR/memdetect.sh"
parseJavaArgs --xmx=2g --xms=32m --mode=fixed "$@"
setEnvironment
if [[ $suite == test.prok.StopTailCorpusTest ]]; then
 exec java $EA $EOOM $SIMD $XMX $XMS -ea -XX:ActiveProcessorCount=1 -cp "$CP" "$suite" "$fixture"
else
 exec java $EA $EOOM $SIMD $XMX $XMS -ea -XX:ActiveProcessorCount=1 -cp "$CP" "$suite" "$@"
fi
