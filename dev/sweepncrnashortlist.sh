#!/bin/bash
# Iterate shortlist policies through the real caller in one JVM.
usage(){
cat <<'HELP'
Usage: dev/sweepncrnashortlist.sh caller=callgenes.config reference=inputs.tsv familymap=rules.tsv policies=policies.tsv out=prefix

Runs whole-input generic-RNA policy iteration with the native caller and reference matcher.
The caller config supplies input FASTAs, enabled families and worker count (t=N).
References use the literal TSV header input_id, fasta, gff. Use reference.tsv
from dev/writencrnareference.sh (reviewed dispositions) or dev/buildncrnareference.sh
mode=tiers (bound CM evidence); see those commands' help for their complete inputs.
Source annotation GFFs are not directly accepted as reference GFFs.
Policies are a TSV: policy, topn, floor for adaptive families, or policy, topn,
minhits for fixed-mode families. Fixed minhits must be integral. Row zero must
match the selected family's packaged baseline; its index mode never changes.
This is an experimental instrument, not a production cutoff selector.

joint=f         With joint=t, vary complete tuples for every member of the
                  selected family. Aggregate cobalamin, 6S and SRP are accepted.
                  The TSV header is exactly:
                  policy, member, mode, topn, floor, minhits, topfrac, qfrac, patience
                  (tab-separated). Mode is fixed or adaptive. Rows for a policy
                  are contiguous; all policies use the same member order, matching
                  the loaded caller order. The first policy preserves every
                  packaged control, including inactive values. Subsequent policies
                  may change mode, fractions and alignment patience. Patience -1
                  retains legacy exits. The first policy must match the loaded
                  packaged patience, including a nonlegacy packaged budget.
                  Nonnegative patience counts eligible non-improving alignments
                  after strict raw-identity improvement; a latched eligible pass
                  enables stopping. Zero stops on the first pass; 2147483647
                  walks the entire primary shortlist. Each primary walk resets;
                  rescue retains its legacy first-pass exit without this budget.
                  Fixed0 and adaptive0 differ for a
                  single model with no shared keys; they are not interchangeable.
                  Joint reports sum work across all selected members, retain the
                  shipped baseline for reference changes, and report other gene
                  counts separately. Equal counts do not prove equal loci.
target=tpp       One registered generic caller, e.g. sam1 or cobalamin_rf00174.
                  Aggregate targets require joint=t; otherwise select a member.
grading=target   Standalone-family grading. Explicit member targets require
                  grading=union and inclusion of every member in familymap.
                  Union counts are not subtype recovery/exclusive accuracy.
cachewindows=128  Retained exact identity windows per worker; 0 disables caching.
cachelength=4096  Largest retained window, in bases.
reusecandidates=f Reuse unchanged pre-path CDS/non-target prok RNA candidates.
                  Every selected member's claims/windows and the global path rerun.
overlap=0.5       Reciprocal scored-reference overlap; masks retain native semantics.
savegff=f         Save complete per-input/per-policy GFFs for bounded verification.
savestages=f      Save .inputN.policyP.stages.tsv with the existing stage schema.
                  The reference binding and policy table name N/P exactly; rows
                  cannot interleave across inputs or policies. Worker-local sinks
                  are rebound before each policy and closed after that input.
                  No shared ncrnadiag/r58diag output may also be configured.
                  Trace files can be large; use only on a bounded diagnostic.
coverage=f        Add root-relative last-loss/first-gain input counts to reports.
members=          Optional frozen family/member/reference_model map, paired with loci=.
loci=             Full qualified support inventory, including unchanged parents.
                  Enables coverage and writes .coverage.tsv matched/scored parent
                  counts for EVERY input/policy/member, including zeroes. Marginal
                  support never confers coverage; member counts are nonadditive.
out=prefix       Fresh .summary.tsv, .inputs.tsv, .changes.tsv and .target.tsv reports.

Cached sweep time is not ordinary caller time. Current alignment counts retain
existing target-scavenger scope. Enable -Dbbtools.diagnostic.ncrnaCpuOut=FILE
in JAVA_TOOL_OPTIONS for actual kernel/model-base, HBM and worker-CPU columns.
Assertions are enabled. Set -Xmx explicitly for the chosen worker count.

Minimal standalone SAM-I example for this version's packaged policy:
  dev/sweepncrnashortlist.sh caller=caller.config reference=ref/reference.tsv \
    familymap=familymap.tsv target=sam1 grading=target joint=t \
    policies=policies.tsv out=trial -Xmx4g
caller.config contains one argument per line (use an absolute original FASTA path):
  in=/absolute/genome.fa
  sam1=t
  t=1
All reference input IDs/FASTA paths must bind that same original input.
familymap.tsv below maps every possible default/opt-in generic caller; include
selects the scoring family, not which callers run. Fields are literal tabs;
family/model_prefix/display_label are nonempty, include is true/t or false/f.
Every enabled generic model must match exactly one model_prefix.
family	model_prefix	display_label	include
RNaseP	rnasep	RNase P	false
tmRNA	tmRNA	tmRNA	false
SRP	srp_small	SRP RNA (small)	false
SRP	srp_large	SRP RNA (large)	false
6S	6S_RF00013	6S	false
6S	6S_RF01685	6S-Flavo	false
FMN	FMN	FMN	false
SAM_I	SAM1	SAM-I	true
TPP	tpp	TPP	false
glycine	glycine	Glycine	false
cobalamin	cobalamin_RF00174	cobalamin_rf00174	false
cobalamin	cobalamin_RF01482	cobalamin_rf01482	false
cobalamin	cobalamin_RF01689	cobalamin_rf01689	false
policies.tsv (literal tabs; first row preserves every packaged field):
policy	member	mode	topn	floor	minhits	topfrac	qfrac	patience
baseline	sam1	fixed	238	0	1	0	0	16
probe	sam1	fixed	238	0	1	0	0	-1
topn>=1; floor/topfrac/qfrac are finite and >=0; minhits is integral and >=0.
Inactive fields remain explicit. The probe is an experimental mapping/walk
comparison, not an emulation of a different complete public caller build.
For coverage, dev/buildncrnareference.sh mode=inventory produces the .loci.tsv.gz
support inventory; use its reviewed family/member/reference_model input table
as members=. Supply the matching witnesses/events and catalog described in
that producer's help. A list of only changed loci is not a complete inventory.
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
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SweepNcrnaShortlist "$@"
