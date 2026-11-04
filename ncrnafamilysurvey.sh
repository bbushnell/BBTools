#!/usr/bin/env bash
usage(){
cat <<'EOF'
NcrnaFamilySurvey: annotated RNA territory and genome presence for bacteria/archaea.
This reads existing annotations; it neither calls genes nor proves biological absence.

Required paths:
map=            Three TAB-separated columns, no header: domain, phylum, basename.
                Domain is bacteria or archaea; basenames must be globally unique.
root=           GFF files are root/domain/basename.gff.gz. Use an absolute root when
                the generated genome roster will feed the annotated-corpus extractor.
outprefix=      Fresh output prefix beneath an existing writable directory.

Selection and output:
expected=39979  Exact complete map count. Set this to YOUR map size; the default is
                a historical corpus size, not an inferred count or universal constant.
limit=0         Zero selects all genomes. A positive limit cycles through sorted
                domain-qualified phyla and names; it is deterministic, not random.
loci=f          Also write individual reconstructed loci when true. Their row order
                is not fixed; exact duplicates merge and retain the earliest source line.
t= / threads=   Defaults to native available-core detection. Set a positive budget
                to fit your allocation. More workers can require more memory.
-Xmx8g -Xms64m Default heap settings; native JVM heap/stack overrides are supported.
                For example: -Xmx4g -Xms64m -Xss2m. Exit-on-OOM is enabled by default.
outstream=     Tool diagnostics (default stderr); JVM/native I/O errors may use stderr.

Outputs: .families.tsv, .phyla.tsv, .genomes.tsv, .summary.tsv, .aliases.tsv,
.excluded_regulatory.tsv; loci=t adds .loci.tsv. Territory sums each locus's interval
union. Presence denominators include selected genomes with zero annotated hits.
Rfam identities, text labels and unknown annotations stay distinct; text labels alone
are not biological family definitions. Fungal/plant/animal inputs are not supported.

config= expands once, in argument order; later settings override earlier settings.
Nested config= directives are rejected. outstream= may appear in this one config.
Help is handled before config expansion. JVM settings stay on the command line.
Native numeric suffixes are supported for counts and threads.
ow=f/app=f are accepted; true and unsupported options are rejected. Existing outputs
are never overwritten. This survey does not read the paired genome FASTA.
EOF
}
if [[ $# -eq 0 || "$1" = -h || "$1" = --help ]];then usage;exit 0;fi
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=fixed -eoom "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaFamilySurvey "$@"
