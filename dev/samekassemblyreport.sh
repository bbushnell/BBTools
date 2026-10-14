#!/bin/bash
# Summarizes collected same-K five-arm QUAST reports; does not tune networks.
# in=collected_directory out=NEW_directory expected=153
# Optional arms=baseline,oracle,nn_three,nn_four,ordinary; exactly five unique names.
# Baseline must be first and ordinary last, preserving paired-comparison semantics.
# Undefined NGA50 is reported as zero with an explicit defined flag/count.
# fold=validation (default) or test must match every manifest row.
# timings=t additionally reads runs/<genome>/{live,ordinary}/{elapsed.tsv,host.txt}.
# timingexclude=g0004 explicitly excludes named trials from timing only.
# Speed summaries require matched hosts and both wall times >=10 seconds;
# timings.tsv retains every pair and its inclusion/exclusion reason.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKAssemblyReport "$@"
