#!/bin/bash
usage(){
cat <<'HELP'
Usage: dev/optimizencrnashortlist.sh caller=callgenes.config reference=inputs.tsv familymap=rules.tsv initial=policies.tsv grid=grid.tsv identity=identity.tsv target=tpp grading=target targetrecall=1 out=search_directory

Experimental local search over coupled shortlist controls. One native sweep runs
each round; every genome keeps normal candidate generation and path selection.
No production cutoff is selected automatically. Resources, seed/window rules,
identity/HBM/length thresholds, boundary and path settings remain fixed. Mapping
and verifier-walk patience vary; their effects on final calls remain measured.

initial=        Complete joint-policy TSV containing exactly two policies:
                the shipped root, then an explicit adaptive starting center.
                Header: policy, member, mode, topn, floor, minhits, topfrac, qfrac, patience
                (tab-separated). Fixed-to-adaptive equivalence must already be
                established by whole-caller testing, not assumed by this tool.
grid=           TSV header: member, knob, values. Supply exactly five rows per
                member: topn, floor, topfrac, qfrac, patience. Values are ascending comma-
                separated choices and must include each starting value. A single
                value freezes an axis. Fine search uses a separate finer grid.
                Patience counts span-eligible alignments since strict raw-identity
                improvement, including before any pass. After a pass is latched,
                stop at this count even when a later attempt fails. Zero stops on
                the first eligible pass; 2147483647 walks the whole shortlist.
                -1 preserves every legacy exit and is outside the search grid.
                The root retains its actual packaged patience, which need not
                be -1. Each primary model walk resets
                its budget; rescue retains its legacy exit and does not use patience.
                Earlier stopping can change selected models, endpoints and path
                competition; paired losses/gains remain in the native reports.
target=         Registered generic family; cobalamin/6S/SRP select whole groups.
grading=        target for standalone families; union for complete groups.
targetrecall=   Required fraction (0,1] of shipped-baseline TP to recover/retain.
                This is relative recovery, not absolute biological recall.
cost=alignments Work objective: alignments uses kernel_logical_alignments;
                modelbases uses aligned_model_bases, summed model length over
                logical alignments. Model bases do not include query length or
                measure DP cells or whole-caller time. Explicit alignments keeps
                the old contract; modelbases binds its unit into the resume contract.
maxlost=-1      Optional maximum count of lost shipped reference IDs. -1 disables
                this constraint; restored TP totals alone do not preserve loci.
maxlastloss=-1  Optional maximum original inputs losing their last matched scored
                target detection relative to shipped. Zero vetoes every such loss.
                First detections gained elsewhere never refund this budget.
maxmemberlastloss=-1 Same limit applied independently to every qualified reference
                member. Member support overlaps; these are not subtype precision.
                Both coverage limits default off and do not change caller output.
members=        Frozen family/member/reference_model map; required with member protection.
loci=           Frozen full reference inventory, including unchanged parents and
                scored_members/marginal_members. Requires member protection.
membersha80=    Declared 20-hex identity of members=. The external recipe verifies it.
locisha80=      Declared 20-hex identity of loci=, authenticated externally, not rehashed
                by this driver. Both support pins are required for member protection.
                The literal small member map and support declarations bind resume.
maxrounds=100   Finite round limit. A local stop is not a global-optimum claim.
heapm=4096      Child JVM heap in MiB; parent heap is set separately with -Xmx.
cachewindows=128  Exact identity windows retained per member per worker.
cachelength=4096  Largest retained window in bases.
reusecandidates=t Reuse unchanged callers' pre-path candidates; all varied members
                  and global path selection execute anew for each policy.
overlap=0.5     Native reciprocal reference-overlap rule, unchanged throughout.
resume=f        Reuse completed rounds only under the same saved contract.
                An incomplete round is reported, never overwritten or rerun.
out=            Fresh search directory (parent exists). No whitespace, quotes or
                backslashes in this path; child CPU profiling uses JVM options.

identity= is a headerless key/value TSV with source_commit (40 lowercase hex),
build_seal, resources_sha80, inputs_sha80 and references_sha80 (20 lowercase hex).
The external run recipe verifies these pins and holds the build lock. This tool
binds the declarations, expanded caller config, manifests, grid and objective;
it does not authenticate runtime or bulk data. Use a controlled JVM environment.

Root-relative lost/gained reference rows and per-input metrics remain in each
round's native reports. Work is summed in the selected unit over all selected
members. Cached CPU is reported separately and is not production runtime.
Below target, loss-free savings/free gains take priority, then recall gained per
extra work unit. At target, keep the target while reducing work cost.
Ties prefer higher TP, lower work cost, fewer lost shipped loci, then fewer
apparent FP, then declared grid order. Lost loci and apparent FP are unbounded
late tie-breakers unless an explicit loss budget applies. Equal outcomes stop.
Results report apparent FP and other gene counts; equal counts do not prove loci
are unchanged. Coarse/fine/protected input selection belongs to the run recipe.
Enabled coverage adds last/first-input counts to the native sweep summary,
coverage_baseline.tsv at the search root, and explicit constraint/target status
plus per-round coverage reasons. Member protection also requires the native
sweep.coverage.tsv with ALL input/policy/member counts, including zeroes and
unchanged matches. Numeric member vectors follow members= order for the graded family.
An infeasible start is reported, never silently reset or relaxed. Resume checks
both the original coverage vector and saved tables, without another caller run.
Omitting both coverage budgets preserves the previous files and result schema.

Input preparation and minimal SAM-I search:
  Use reference.tsv from dev/writencrnareference.sh (reviewed dispositions) or
  dev/buildncrnareference.sh mode=tiers. Their help gives the complete producer inputs.
  The reference manifest header is input_id, fasta, gff (literal tabs and absolute
  paths). Reference GFFs need the producer's explicit input/status/evidence fields.
  caller.config has one flag=value per line, for example:
    in=/absolute/genome.fa
    sam1=t
    t=1
  familymap.tsv header is family, model_prefix, display_label, include, with
  literal tabs. include is true/t or false/f. Map every enabled generic model
  exactly once, including explicitly excluded families; set only SAM_I included.
  dev/sweepncrnashortlist.sh --help prints a complete ready-to-copy SAM-I family map
  for the default and opt-in prokaryotic callers, plus the identical input rules.
  For optional member protection, dev/buildncrnareference.sh mode=inventory produces
  .loci.tsv.gz; members= is its reviewed family/member/reference_model input.
  Use the same full reference/catalog/witness/event set, not a changed-locus list.

initial.tsv (literal tabs; root is the packaged point, start explicitly adaptive):
policy	member	mode	topn	floor	minhits	topfrac	qfrac	patience
baseline	sam1	fixed	238	0	1	0	0	16
start	sam1	adaptive	100	1	1	0.1875	0	16
grid.tsv (literal tabs; the one-value axes intentionally stay fixed):
member	knob	values
sam1	topn	80,100,128
sam1	floor	1
sam1	topfrac	0.1875
sam1	qfrac	0
sam1	patience	8,16,32
Then, with the declared identity.tsv described above:
  dev/optimizencrnashortlist.sh caller=caller.config reference=ref/reference.tsv \
    familymap=familymap.tsv initial=initial.tsv grid=grid.tsv identity=identity.tsv \
    target=sam1 grading=target targetrecall=1 maxlost=0 out=new_search -Xmx4g
This bounded example illustrates file construction, not a recommended cutoff.
HELP
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.OptimizeNcrnaShortlist "$@"
