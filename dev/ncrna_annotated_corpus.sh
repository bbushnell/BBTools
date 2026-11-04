#!/usr/bin/env bash
# Extract an Rfam population using the survey's source-locus definition.
usage(){
cat <<'EOF'
Extract a complete annotated Rfam population from a frozen survey roster.
Usage: ncrna_annotated_corpus.sh roster=survey.genomes.tsv family=RF00174
       expectedgenomes=39979 outprefix=corpus t=8

roster=           Existing thirteen-column native survey genome table, with absolute
                  paired .gff.gz paths; corresponding genomes end in .fna.gz.
family=           Exact RF followed by five digits; source flags are not independent truth.
expectedgenomes=  Required complete roster size, checked even when limit selects a pilot.
outprefix=        Fresh prefix for .all.fa.gz, .full.fa.gz, .contexts.fa.gz,
                  .loci.tsv.gz, .genomes.tsv and .summary.tsv. Never overwrites.
expectedloci=-1   Optional exact extracted raw-locus count; -1 disables this check.
limit=0           First N genomes for an explicitly labeled pilot; 0 retains the roster.
flank=200         Genomic context on either side, clipped at contig edges; range 0..2000.
t=8               Genome workers. JVM heap flags such as -Xmx8g use native launcher setup.
gccontig=f        Add contig_gc only to context FASTA headers when true. Uses the entire
                  original contig before cropping or strand reversal: G/C over defined
                  bases, U as T, ambiguous bases ignored; zero if no base is defined.
                  Float encoding matches inference. It is not locus or window GC.
                  Default false leaves all existing output formats unchanged.

Exact genes concatenate reconstructed intervals in biological orientation; contexts
remain contiguous and retain real flanks. All raw loci remain in all/context/ledger
outputs; full.fa contains only simple annotation-full copies. No source downloads,
resource training or scheduler submissions are performed by this extractor.
EOF
}
if [[ $# = 0 || "$1" = -h || "$1" = --help ]];then
  usage
  exit 0
fi
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=fixed -eoom "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaAnnotatedCorpus "$@"
