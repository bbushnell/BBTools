#!/usr/bin/env bash
usage(){
cat <<'HELP'
Written by Keqing.

Calibrates a protein-family classifier using complete negative panels.
Loads the classifier once, scores with one worker pool per round, and shares its
unchanged index and features across rounds. Optional preparation builds the initial
classifier in memory from pinned member inputs. Positive sequences never enter calibration.

Usage: familycalibrate.sh in=classifier.config neg=k5.faa.gz,shadow.faa.gz \
         negnames=k5,shadow targetfpr=0.0001 rounds=4 out=calibrated t=32
   or: familycalibrate.sh prepare=preparation.config out=calibrated t=32

in=              Existing classifier recipe, or directory containing classifier.config.
                 Retains that package's stored matrix and scoring policy.
                 Compact SHORT/TAILS family files are not supported for recalibration.
prepare=         Pinned family_preparation_v1 manifest, instead of in/neg/negnames/negsha80.
                 Supplies members, fixed consensuses, labeled negative panels, explicit
                 geometry, retention, memory/output bounds and deployment settings.
                 Passes references and modal profiles directly to calibration in memory.
                 Writes final references/config once after the last evaluated round,
                 retaining profiles/, raw_calibration.tsv.gz, preparation.tsv and
                 preparation_inputs.tsv beside the final package; no initial/ copy.
                 Explicit alignmentmode/alignmentminspan bind the measured geometry.
                 New-metric preparation also requires heuristic=hitfourth_h020_c150,
                 pinned lengthmedians/familyhydro and minqueryaa/maxqueryaa in the
                 manifest. Unavailable alignment bounds remain explicit in the audit.
                 New preparation defaults to POSSUM55 from resources/bbgap/POSSUM55.tsv;
                 set any explicit matrix selector and its resource/binding in the manifest.
neg=             Comma-separated complete negative FASTA/FASTQ panels.
negnames=        Optional distinct panel labels; positive is reserved.
negsha80=        Optional input identities (20 hex digits), in panel order.
                 Structural mode records supplied identities as inherited and leaves
                 omitted identities unmeasured. Explicit SHA80 mode checks supplied
                 fingerprints and measures missing ones.
targetfpr=0.0001 Maximum accepted fraction, inclusive; finite value from 0 to 1.
targetaccepts=   Optional exact count instead of targetfpr; requires perclass=t,
                 one existing-package negative panel and expectedpanelqueries.
marginscale=1    Positive finite evaluation margin for existing-package repair.
                 Raises effective gates at this margin, then converts back to stored
                 upper gates. Use the same margin when evaluating the exported package.
verification=structural  Default: validate input structure, counts and stability
                 without computing content hashes. Ordinary existing-package panels
                 may omit identities and expected counts. Use verification=sha80
                 for exact content-hash checking; inherited labels are not byte proof.
querygc=         Optional absolute ordered query-GC manifest; mandatory for a GC-ranked
                 package. Requires structural existing-package calibration, without
                 prepare or combined heldout scoring. TSV columns: panel, path, queries,
                 source_kind, query_gc, annotation_rows, sidecar, contigs, genome_manifest,
                 genomes. One row per negative panel; kinds annotation/genomic_shadow/none.
                 Unused paths are '-' and unused counts 0. Load once and share across
                 rounds; preserve all panels and explicit neutral GC. Writes final
                 gc_coverage.tsv and per-round gc_rounds.tsv, without per-query files.
perclass=f       Apply target to pooled counts. With t, every class must meet it.
rounds=4         Maximum evaluations, including the initial parent; range 1..32.
expectedinitialaccepts= Optional comma-separated accepted counts in negative-panel
                 order. A round-one mismatch fails before any gate update.
expectedpanelqueries= Optional exact count for every panel/round, or a comma-separated
                 vector of exact counts in panel order. Query-GC and exact-count
                 repair require declared populations. Without this assertion,
                 calibration compares later observed counts with the first round.
out=             New output directory; partial and completed runs are never overwritten.
t=               Scoring workers; threadsin=1 sets input workers per panel.

Each round observes one immutable parent across all panels. Identity winners
raise their family gate to score+1 percentage point; BLOSUM winners use score+0.005.
Maxima merge before one upward-float-rounded update. Lower gates and ranking stay fixed.

Outputs: classifier.config, references.tsv.gz for prepare= or changed gates,
summary.tsv, family_false_accepts.tsv.gz, rounds.tsv, cohorts.tsv, calibration.tsv,
inputs.tsv, interpretation.txt and PASS. FROZEN exists only for convergence.
PASS means the calculation completed, not that it converged: inspect status for
converged, round_limit, ddl_floor_exceeds_target or no_progress. Nonconvergence
retains the last evaluated package. These are calibration errors, not cold FPR.
Shared pinned ranking files are referenced rather than copied into each output.

dependencyjar= accepts colon-separated existing dependency jars, after modulejar
and before the runtime on the classpath. Default classpath is unchanged.
classpathfile= reads one literal classpath line, preserving provider order.
Entries must be existing absolute files or directories. Cannot combine with
modulejar= or dependencyjar=. library= still selects JVM setup scripts.
Standard config= and -Xmx options are supported.

Development-module-only validation selectors take one fresh absolute fixture directory:
gcpenaltytest, gcinputstest, gcworkertest, gcpackagetest, querygctest,
gcpaneltest and gcdrivertest. These synthetic fixtures require a separately supplied
test module on the classpath; their test classes are not included in the shipped
selection. Supply the module with modulejar= or include it in classpathfile=.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]];then usage;exit 0;fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
LIBRARY=$DIR;MODULE=;DEPENDENCY=;ARGS=();JVM_ARGS=()
CLASSPATH_FILE=;EXACT_CP=false;MODULE_SET=false;DEPENDENCY_SET=false
for arg in "$@";do
 case ${arg,,} in
  library=*) LIBRARY=${arg#*=};;
  modulejar=*) MODULE=${arg#*=};MODULE_SET=true;;
  dependencyjar=*) DEPENDENCY=${arg#*=};DEPENDENCY_SET=true;;
  classpathfile=*)
   [[ $EXACT_CP == false ]] || { printf 'classpathfile may be supplied only once\n' >&2;exit 2; }
   CLASSPATH_FILE=${arg#*=};EXACT_CP=true;;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da) printf '%s\n' 'Assertions must remain enabled' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done
set -- "${ARGS[@]}"
DIR=$LIBRARY
if [[ $EXACT_CP == true ]];then
 [[ $MODULE_SET == false && $DEPENDENCY_SET == false ]] || {
  printf 'classpathfile cannot combine with modulejar or dependencyjar\n' >&2;exit 2;
 }
 [[ -n $CLASSPATH_FILE && -f $CLASSPATH_FILE && -s $CLASSPATH_FILE ]] || {
  printf 'classpathfile requires an existing nonempty file\n' >&2;exit 2;
 }
 mapfile -t cp_lines < "$CLASSPATH_FILE" || exit 2
 [[ ${#cp_lines[@]} == 1 ]] || { printf 'classpathfile requires exactly one line\n' >&2;exit 2; }
 CP=${cp_lines[0]}
 [[ -n $CP && $CP != :* && $CP != *: && $CP != *::* && $CP != *$'\r'* ]] || {
  printf 'classpathfile has an empty provider or invalid line ending\n' >&2;exit 2;
 }
 IFS=: read -r -a cp_entries <<< "$CP"
 for path in "${cp_entries[@]}";do
  [[ $path == /* && ( -f $path || -d $path ) ]] || {
   printf 'classpathfile provider must be an existing absolute file or directory: %s\n' "$path" >&2;exit 2;
  }
 done
else
 if [[ -f $DIR/bbtools.jar ]];then CP=$DIR/bbtools.jar;else CP=$DIR/current;fi
 if [[ -n $DEPENDENCY ]];then
  [[ $DEPENDENCY != :* && $DEPENDENCY != *: && $DEPENDENCY != *::* ]] || {
   printf 'dependencyjar must not contain an empty classpath provider\n' >&2;exit 2;
  }
  IFS=: read -r -a dependencies <<< "$DEPENDENCY"
  for path in "${dependencies[@]}";do [[ -f $path ]] || exit 1;done
  CP=$DEPENDENCY:$CP
 fi
 if [[ -n $MODULE ]];then [[ -f $MODULE ]] || exit 1;CP=$MODULE:$CP;fi
fi
CLASS=prot.FamilyCalibrate
if [[ ${1:-} == test ]];then CLASS=test.prot.FamilyCalibrateTest;shift;fi
if [[ ${1:-} == gcdrivertest ]];then CLASS=test.prot.FamilyGcDriverTest;shift;fi
if [[ ${1:-} == gcpenaltytest ]];then CLASS=test.prot.FamilyGcPenaltyTest;shift;fi
if [[ ${1:-} == gcinputstest ]];then CLASS=test.prot.FamilyGcInputsTest;shift;fi
if [[ ${1:-} == gcworkertest ]];then CLASS=test.prot.FamilyGcWorkerTest;shift;fi
if [[ ${1:-} == gcpackagetest ]];then CLASS=test.prot.FamilyGcPackageTest;shift;fi
if [[ ${1:-} == querygctest ]];then CLASS=test.prot.FamilyQueryGcTest;shift;fi
if [[ ${1:-} == gcpaneltest ]];then CLASS=test.prot.FamilyGcPanelTest;shift;fi
if [[ ${1:-} == preparationtest ]];then CLASS=test.prot.FamilyPreparationTest;shift;fi
if [[ ${1:-} == metricpreparationtest ]];then CLASS=test.prot.FamilyPreparationMetricTest;shift;fi
if [[ ${1:-} == pilotinput ]];then CLASS=test.prot.FamilyPreparationPilotInput;shift;fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=auto "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" "$CLASS" "$@"
