#!/usr/bin/env bash
usage(){
cat <<'HELP'
Written by Keqing.

Evaluates positive and negative protein panels using one loaded family classifier.
Uses the existing classifier scoring path; does not rewrite gates or change families.

Usage: familyeval.sh in=classifier.config pos=heldout.faa.gz \
         neg=k5.faa.gz,shadow.faa.gz negnames=k5,shadow out=results t=32

in=          Classifier package recipe, or directory containing classifier.config.
             The recipe contains existing FamilyTailCascade callable key=value
             settings, including familyreferences and all scoring policy choices.
             Component paths may be relative to the recipe. Optional sha80 labels
             are inherited in structural mode; exact SHA80 mode checks supplied
             component pins and measures missing component pins.
pos=         Complete unpaired positive FASTA/FASTQ panel.
neg=         Comma-separated complete negative panels. Each has its own denominator.
             Supply either or both. Negatives alone can select a policy before positives are scored.
negnames=    Optional unique labels, one per negative panel. Default negative_1, etc.
possha80=    Optional expected positive-file fingerprint (final 20 hex digits).
negsha80=    Optional comma-separated negative-file fingerprints, in panel order.
verification=structural  Default: validate input structure, counts and stability
             without computing content hashes. Supplied identities are inherited;
             omitted identities remain unmeasured. Use verification=sha80 for
             exact content-hash checking; inherited labels are not byte proof.
poscount=    Optional expected positive-panel population in structural mode.
negcounts=   Optional expected populations in negative-panel order. Ordinary panels
             may omit these counts; scoring reports the observed denominators.
             Query-GC requires exact counts. Omitted expected counts are reported
             as NA in input_populations.tsv, not as measured populations.
evaluationgc=package  Default: retain the package's GC ranking and filtering.
             Set evaluationgc=off explicitly to omit both for this evaluation.
             Does not alter the package, its header or caller GC policy. Records
             the override in evaluation_policy.tsv; the emitted reference/parent
             recipe still identifies the original package. Ordinary evaluation
             only; cannot combine off with querygc, paired ranking or saved floors.
querygc=     Optional absolute ordered query-GC manifest. Requires verification=structural
             and exact poscount/negcounts. Mandatory for a GC-ranked package
             unless evaluationgc=off is explicitly selected; absence never disables GC.
             TSV columns: panel, path, queries, source_kind, query_gc, annotation_rows,
             sidecar, contigs, genome_manifest, genomes. One row per existing panel;
             source_kind is annotation, genomic_shadow or none. Unused paths are '-'
             and unused counts are 0. Missing measured GC is explicitly neutral;
             no panel is dropped. Resident bindings load once, including benchmark
             repeats. Writes gc_coverage.tsv; no per-query intermediate is written.
floors=      Optional frozen fitted-family table (raw TSV or gzip/BGZF).
floorssha80= Required table label: lowercase 20 hex digits, or literal unmeasured
             in structural mode. Also require insha80= for the original
             classifier.config, fitting= for its fitting-cohort manifest,
             fittingsha80= and explicit marginscale= matching the saved table.
             Parent/cohort labels must match the saved table metadata exactly;
             unmeasured is not a wildcard. Every panel requires a declared label.
             Structural mode checks observable file stability, not byte identity.
             verification=sha80 requires 20-hex labels matching the actual bytes.
             This mode loads one bound classifier, permits one complete evaluation,
             and rejects ranking/cap/admission overrides, paired mode and repeats.
             Query bounds are fixed at40..99. It never refits floors/support.
             The output parent_classifier.config is only the parent recipe;
             evaluation_policy.tsv and inputs.tsv retain floor/cohort bindings.
             No classifier.config is emitted that could silently omit the floors.
out=         New result directory; existing directories are never overwritten.
t=           Scoring workers. All panels share the same loaded classifier and pool.
threadsin=1  Input workers per panel.
top=         Optional positive shortlist-cap override; parent recipe remains unchanged.
ddlreject=   Optional fixed raw-hit admission override, e.g.6,5,4. Copied into
             the loaded policy and resolved output recipe; parent bytes stay unchanged.
             May be combined with heuristic ranking. Length policies require anchor6.
admission=fixed  Keep the package's fixed raw-hit floor, or select family_mean:
             4 hits for family mean<65aa,5 for65<=mean<100,6 formean>=100.
             Uses exact retained-member count/residue totals already in the reference;
             no rounded mean, query length, consensus length or new export is used.
             Requires the new heuristic, top50/P4/O2, exhaustive indexing and
             disabled shortcut. Existing family gates and cosine arrays stay fixed.
             Legacy short_scaled remains a separate explicit comparison:
             ceil(6*(L-4)/95) for normalized40..99-aa exact-K5 queries.
             This experiment requires heuristic=hitfourth_h020_c150, top50/P4/O2,
             exhaustive indexing and disabled sketch shortcut. Other lengths
             fail; no records are filtered. Family gates and cosine arrays stay fixed.
admissionstats=f  Write admission_summary.tsv with raw admitted-pool counts
             overall and in requested lengthbins. Enabled by both length policies;
             use admissionstats=t with fixed admission for the matched control.
             floor_min/floor_max report configured family-floor bounds for family_mean.
stagecpu=f  Write stage_cpu.tsv with opt-in worker CPU counters split into
             sketch/index lookup, shortlist/ranking, alignment/acceptance and other work.
             This is diagnostic timing only and does not change scoring.
heuristic=off  Use legacy ranking, or hitfourth_h020_c150 for
             hits*fourth_root(length ratio)+1.5*cosine+0.2*member-KD similarity.
             The heuristic requires pinned lengthmedians/familyhydro tables and
             explicit minqueryaa/maxqueryaa. It preserves existing cosine arrays.
familyhydro= Absolute retained-member KD table; familyhydrosha80= is required.
minqueryaa=  Inclusive normalized query-length bounds, supplied together with
maxqueryaa=  maxqueryaa, within1..200. Validate only: no filtering, cropping or
             sampling. One conventional terminal '*' is normalized as before.
             Use1/200 for the complete capped-tail cohort, not200/200.
lengthbins=  Optional increasing positive integer boundaries, e.g.40,60,80,100.
             Writes length_summary.tsv from the same scoring pass, using normalized
             peptide lengths and [lower,upper) bins. Includes underflow/overflow
             and empty bins; NA upper bound means unbounded, NA rates mean no queries.
             No filtered FASTAs or additional classifier loads are required.
             Reporting only; ordinary single-policy evaluation, not paired mode.
marginscale=1  Evaluation-only identity/BLOSUM acceptance margin scale; finite and >=0.
             For each metric, base=max(member lower,global floor),
             current=max(frozen upper,global floor), threshold=base+scale*(current-base).
             One preserves the original threshold exactly; zero uses the base.
             Values below one loosen positive margins; above one tighten them.
             Strict greater-than comparisons use double thresholds. Zero margins stay fixed.
             Member rejection, ranking, scalar gates and DDL/cosine shortcuts do not change.
             Frozen reference files remain unchanged; the scale is recorded separately.

paired=f     Opt-in joint OFF/multiply/add1 length-ranking comparison.
             Paired mode requires positive-only OR negative-only panels, never both.
lengthmodel= Absolute frozen epoch-model path; lengthmodelsha80= is required.
lengthmedians= Absolute family-median table for paired or heuristic ranking;
             lengthmedianssha80= is required. The two ranking modes cannot mix.
pairedscales= Three explicit absolute margin scales in OFF,multiply,add1 order.
             Use instead of marginscale. Select scales from negatives in separate
             invocations and freeze before positives; this command never selects them.
             Shared retrieval and union alignment support paired_transitions.tsv and
             oracle_membership.tsv. Oracle bits are OFF=1,multiply=2,add1=4;
             zero means an empty union. The oracle is restricted to the union of
             strict pre-scalar top50 lists, not exhaustive family truth.

Outputs summary.tsv (counts, TPR/FNR, FPR and false accepts per million per class),
family_false_accepts.tsv.gz (terminal-winner counts for every live family, including
zero-count rows), resolved classifier.config, evaluation_policy.tsv, inputs.tsv, interpretation.txt and
timing.tsv. Empty reference models are not callable and do not enter family rows.
Ordinary mode also writes scoring_time.tsv: wall_seconds, process_cpu_seconds,
queries and worker_threads. This interval includes streaming, worker setup,
classification and reduction; excludes model/index loading, pin verification
and output writing. It is not whole-process CPU or pure alignment time.
TPR means acceptance by any family, not correct-family assignment. Reusing a
calibration negative panel does not make its error rate independent validation.
Empty panels, invalid package settings and nonfinite metrics fail with a reason.

library= selects the runtime; modulejar= overlays the evaluator.
dependencyjar= accepts colon-separated existing dependency jars, after the module
and before the runtime on the classpath.
classpathfile= reads one nonempty literal classpath line, preserving provider order.
Every colon-separated entry must be an existing absolute file or directory.
Cannot combine with modulejar= or dependencyjar=; library= still selects JVM
setup scripts. No runtime provider is appended when classpathfile= is supplied.
Standard config= options and -Xmx memory control are supported.

Native structural fixture: familyeval.sh structuraltest /absolute/new-directory
Runs synthetic checks in a fresh directory. Takes that path plus launcher/JVM
options only; the fixture selects its own product settings and worker counts.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then usage; exit 0; fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
LIBRARY=$DIR; MODULE=; DEPENDENCY=; ARGS=(); JVM_ARGS=()
CLASSPATH_FILE=; EXACT_CP=false; MODULE_SET=false; DEPENDENCY_SET=false
for arg in "$@"; do
	case ${arg,,} in
		library=*) LIBRARY=${arg#*=};;
		modulejar=*) MODULE=${arg#*=}; MODULE_SET=true;;
		dependencyjar=*) DEPENDENCY=${arg#*=}; DEPENDENCY_SET=true;;
		classpathfile=*)
			[[ $EXACT_CP == false ]] || { printf 'classpathfile may be supplied only once\n' >&2;exit 2; }
			CLASSPATH_FILE=${arg#*=}; EXACT_CP=true;;
		-xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
		-da) printf '%s\n' 'Assertions must remain enabled' >&2; exit 2;;
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
 if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
 if [[ -n $DEPENDENCY ]]; then
  [[ $DEPENDENCY != :* && $DEPENDENCY != *: && $DEPENDENCY != *::* ]] || {
   printf 'dependencyjar must not contain an empty classpath provider\n' >&2;exit 2;
  }
  IFS=: read -r -a dependencies <<< "$DEPENDENCY"
  for path in "${dependencies[@]}"; do [[ -f $path ]] || exit 1; done
  CP=$DEPENDENCY:$CP
 fi
 if [[ -n $MODULE ]]; then [[ -f $MODULE ]] || exit 1; CP=$MODULE:$CP; fi
fi
CLASS=prot.FamilyEval
case ${1:-} in
	test) CLASS=test.prot.FamilyEvalTest; shift;;
	pairedcountstest) CLASS=test.prot.FamilyPairedCountsTest; shift;;
	pairedtest) CLASS=test.prot.FamilyPairedEvalTest; shift;;
	structuraltest) CLASS=test.prot.FamilyEvalStructuralTest; shift;;
	cascadetest) CLASS=test.prot.FamilyTailCascadeTest; shift;;
	positionaltest) CLASS=test.prot.FamilyTailCascadePositionalTest; shift;;
	shortadmissiontest) CLASS=test.prot.FamilyShortAdmissionTest; shift;;
	floortest) CLASS=test.prot.FamilyHitFloorEvalTest; shift;;
esac
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=auto "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" "$CLASS" "$@"
