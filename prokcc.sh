#!/usr/bin/env bash
# Native assembly-batch completeness/contamination client.
set -o pipefail
usage(){
cat <<'USAGE'
Description:  Estimates completeness and contamination of prokaryotic assemblies.
              One FASTA file is one bin. Models load once for the whole batch.

Usage:  prokcc.sh in=<assembly.fa> out=<quality.tsv> <other arguments>
or
prokcc.sh <assembly.fa> <other arguments>

Examples:
prokcc.sh ecoli.fa
prokcc.sh in=bin.fa out=quality.tsv t=8 -Xmx8g
prokcc.sh in=bin1.fa,bin2.fa out=quality.tsv
prokcc.sh in=/directory/of/bins out=quality.tsv

File parameters:
in=<file>             Assembly FASTA, comma-separated files, or a directory.
                      A bare existing input file is also accepted. Directories
                      include FASTAs in sorted order, without recursion.
out=<file>            Optional TSV data output; omitted by default.
                      Use out=stdout for TSV on stdout; out=null discards TSV.
                      Results are published only after every bin succeeds.
ow=t                  Overwrite existing output. Set ow=f to protect existing files.
config=<file>         Release configuration. Default: resources/prokcc/release.config
                      beside this installation. Relative resource paths use the
                      config directory; in/out paths use the working directory.

Processing parameters:
t=<integer>           Bin workers; defaults to available threads, capped by inputs.
taxaddress=refseq     QuickClade server for taxonomy. Input-header taxonomy is ignored.
normalsearch=t        Uses the normal candidate pool and reranking for taxonomy.
                      Set false for the legacy one-hit search. Normal server
                      search requires a server that acknowledges this mode.
taxdomain=<name>      Bacteria or Archaea; bypasses QuickClade when supplied.
taxphylum=<name>      Optional phylum with taxdomain=. Rows are marked user-supplied.
pgmmode=taxonomy      Gene-caller model selection: taxonomy or default.
passes=1              Gene-calling passes.
deterministic=t        Preserves the reference inference arithmetic.
loadmode=parallel     Resource loading: parallel or serial. Independent of t=.
compositemode=locked  Composite sharing: locked or worker (private copies).
subnetmode=locked     Subnet sharing: locked or worker (private copies).
timings=f             Reports process-phase wall times to stderr when enabled.
verbose=f             Prints resource-loading and inference diagnostics to stderr.
swapnl=f              N is count and L is length by default. Set true to reverse
                      the N/L labels, as in stats.sh. Alias: swapln.
comperrormultiplier=1.0    Positive finite multiplier for completeness error.
contamerrormultiplier=1.0  Positive finite multiplier for contamination error.
selftest=f            Runs small CLI/input fixtures instead of scoring assemblies.

Java parameters:
-Xmx8g                Maximum heap. Allow room for models and per-worker state.
-ea                   Assertions enabled by default.

Model installation:
Download prokcc_v1.2.1.tar from:
https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.2.1.tar
Extract its contents into resources/ to create resources/prokcc/.
The archive contains the matching configuration, models, and tables. Downloads
are manual. See resources/prokcc/README.md for the model and release details.

Output:
Aligned human-readable reports always go to stderr, one block per bin. Without
out=, stdout is empty. Completeness, contamination, their errors, ANI, and coding
density are percentages in the screen reports. GC is a fraction.
The TSV contains completeness/contamination, six raw model heads, predicted
absolute errors, taxonomy provenance, input sha80, reference name/TaxID/ANI,
contig Nx/Lx, genome size, GC, CDS/RNA counts, coding density, and MIMAG tier.
Scores, errors, ANI, GC, and coding density are fractions, not percentages;
error estimates are not confidence intervals. Missing reference metrics are NA.
The final column is per-bin worker time, excluding shared setup and publication.
USAGE
}
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then usage; exit 0; fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")" && pwd) || exit 1
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current/; fi
source "$DIR/javasetup.sh" || exit 1
CLASS=prok.ProkCC
ARGS=(); JVM_ARGS=(); SELFTEST=''; CONFIG_GIVEN=false
for arg in "$@"; do
  case ${arg,,} in
    config=*) CONFIG_GIVEN=true; ARGS+=("$arg");;
    selftest=*)
      [[ -z $SELFTEST ]] || { echo 'Duplicate selftest flag' >&2; exit 2; }
      SELFTEST=${arg#*=}
      case ${SELFTEST,,} in
        t|true) CLASS=prot.MagQCAssemblyBatchTest;;
        f|false) ;;
        *) echo 'selftest must be t or f' >&2; exit 2;;
      esac;;
    --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
    -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
    simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
if [[ $CLASS == prok.ProkCC && $CONFIG_GIVEN == false ]]; then
  DEFAULT_CONFIG=$DIR/resources/prokcc/release.config
  if [[ ! -f $DEFAULT_CONFIG ]]; then
    echo "Missing default model config: $DEFAULT_CONFIG" >&2
    echo 'Download prokcc_v1.2.1.tar from https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.2.1.tar' >&2
    echo "Extract its contents into $DIR/resources/ to create resources/prokcc/, then try again." >&2
    exit 1
  fi
  ARGS=("config=$DEFAULT_CONFIG" "${ARGS[@]}")
fi
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}" || exit 1
setEnvironment || exit 1
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
