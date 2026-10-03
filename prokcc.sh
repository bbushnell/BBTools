#!/usr/bin/env bash
# Native assembly-batch completeness/contamination client.
set -o pipefail
usage(){
cat <<'USAGE'
Usage: prokcc.sh in=bin.fa out=quality.tsv config=release.config t=8 -Xmx8g
       prokcc.sh in=bin1.fa,bin2.fa out=stdout config=release.config
       prokcc.sh in=/directory/of/bins out=quality.tsv config=release.config

Without config=, use resources/prokcc/release.config beside this installation.
That file names the preserved current model; see resources/prokcc/README.md.
Download prokcc_v1.tar from:
https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.tar
Extract its contents into resources/ to create resources/prokcc/. The archive contains the matching
release config, models and tables; the tool does not download them automatically.

One FASTA is one bin. Directory input includes FASTA files in sorted order,
without recursion; comma-separated input retains its order. Files must remain
unchanged during the run. bin_id is the normalized absolute input path.
Networks and assignment resources load once; workers process bins concurrently.
The complete report is published only when every bin succeeds. Output files
must be fresh. t/threads controls bin workers; each native FASTA reader uses
one input thread. Size heap for the resources plus worker-local inference state.

The release config supplies the real six-output model and frozen resources:
net/netsha80, bundle/bundlesha80, familylist/familylistsha80,
subnetmanifest/subnetmanifestsha80, expectedcopytable/expectedcopytablesha80,
subnetpopulations/subnetpopulationssha80, profile/profilesha80,
roster, ref, rolemanifest, core, coveringsets, sidecar, hbmbundle, hbmprovenance.
Relative RESOURCE paths resolve against that config's directory, not the
working directory. in/out paths still resolve against the working directory.
Use one config file. Duplicate/unknown options and dummy models are rejected.

taxaddress=refseq uses QuickClade. taxdomain=Bacteria|Archaea with optional
taxphylum=NAME bypasses the server and marks each row as user-supplied taxonomy.
Input-header taxonomy is never used. A valid no-hit is retained as unknown;
server failures are fatal. pgmmode=taxonomy and passes=1 are defaults.
Assignment uses the release policy BOUNDED_LOOKAHEAD/lookahead=4.
deterministic=t preserves the existing reference inference arithmetic.
compositemode=worker|locked and subnetmode=worker|locked select inference sharing
(both default to locked). Locked mode uses one lazy composite or one
locked dense inference instance per subnet. Formatter inputs and returned heads
remain private. Use worker for independent network copies per bin worker.
loadmode=serial|parallel selects sequential resource loading or three concurrent
loads (HBMs, subnet bundle, composite). Default parallel; bin workers still use t=.
timings=t reports optional phase wall times to stderr; worker phases are summed
across bins. A single-bin run gives an additive split including setup/I/O.

Output retains all six raw model heads and adds two error estimates:
raw error times comperrormultiplier/contamerrormultiplier (default1.0).
Factors must be positive and finite; products are not clipped. Config metadata
errorfitset=UNCALIBRATED, errorfitdate=NA, errorcoverage=NA remain explicit until
a reviewed calibration supplies them. Errors are fraction units, not confidence
intervals. Each row includes taxonomy provenance and the input sha80.
The v3 TSV also includes reference name/TaxID/ANI, whole-record contig Nx/Lx,
genome size, GC/ACGT, CDS/RNA counts, summed coding bp/genome size, and the
RNA-aware extended MIMAG tier. TSV ANI, GC and coding density are fractions.
Missing reference metrics or GC without ACGT are NA. Overlapping CDS can make
coding density exceed 1. Single-assembly input also prints an aligned summary
to stderr, with ANI/completeness/contamination/coding density in percent.
Each row ends with bin_worker_wall_seconds: that bin's taxonomy and calling/
inference worker intervals, including input I/O. It excludes shared setup,
dispatch queues, phase barriers and final report publication. It is always
measured, independently of the optional timings=t process-phase diagnostics.

selftest=t runs only small CLI/input-contract fixtures, not biological validation.
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
    echo 'Download prokcc_v1.tar from https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.tar' >&2
    echo "Extract its contents into $DIR/resources/ to create resources/prokcc/, then try again." >&2
    exit 1
  fi
  ARGS=("config=$DEFAULT_CONFIG" "${ARGS[@]}")
fi
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}" || exit 1
setEnvironment || exit 1
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
