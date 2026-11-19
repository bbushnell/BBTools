#!/usr/bin/env bash
# One exported positional network and one held-out matrix/end; scores saved once.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Score one exported ncRNA positional network on saved held-out jittered trials.
Usage: evaluatencrnapositionalnet.sh in=evaluation.tsv rows=evaluation.rows.tsv \
 net=endpoint.bbnet end=5prime out=summary.tsv

Requires a plain #dims TAB24|69 TAB21 TAB0 matrix and its eleven-column evaluation
sidecar. end=5prime or end=3prime selects the saved sense-coordinate label.
Reports trial counts, exact/within1/within3, mean/RMS endpoint error, and movement
of initially correct ends (percentage denominator: only initially correct trials).
No unique-genome recall, joint-end validity or whole-caller detection is measured.
Raw features pass through serialized network normalization exactly once. Tied
outputs keep the first maximum, matching the caller at zero confidence margin.
The original out= summary retains its format. Two additional outputs are default:
 raw=<out>.scores.tsv             All 21 scores as lossless signed FP32 bit words,
                                 with exact identity rows and a COMPLETE footer.
 comparison=<out>.comparison.tsv  Margins 0 and 0.1, plus doing nothing, including
                                 ungated corrections kept/lost by the margin.
All policies use the same forward pass. The shared caller confidence rule holds
at margin 0.1 unless best beats both current and runner-up by more than 0.1.
Paths are recorded in the raw header; bind input/network fingerprints and runtime
flags in the run manifest. Completion requires a successful process and closed
outputs, not the footer alone. Summaries follow complete input and raw closure.
Append is unsupported. Use ow=t to replace outputs. This tool neither trains nor
selects a winning network.

fullerrors=f   Optional69-input/22-output path with the explicit out-of-range head.
                Use fullerrors=t only with its exact69/22header and matched
                evaluation sidecar. Both margins0/0.1 hold unless the best real
                offset strictly beats current, runner-up AND out-of-range by
                more than the margin. Ties hold even at margin0. Out-of-range
                never means a +11-base move. out= is one full policy table with
                do-nothing beside both margins; comparison= is unsupported here.
                RMS uses ALL actual sidecar errors including distant tails.
                Counts/percentages include exact/within1/within3, correct ends
                moved, reachable wrong ends fixed, unreachable ends held/moved,
                and separate reachable OOR-head wins versus additional OOR vetoes.
                raw= retains all22FP32 score words plus exact signed correction
                and actual policy selections, with a COMPLETE footer. Default
                21-output behavior and output formats remain unchanged.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.EvaluateNcrnaPositionalNet "$@"
