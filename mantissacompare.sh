#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 2026

Description:  Single-bucket sub-NLZ state simulator.  Simulates many independent
trials of a single register (bucket) being fed random hashes, tracking how the
sub-NLZ state bits (history, mantissa, luck, andtissa, nlz2, or combinations)
evolve as cardinality grows.  For each NLZ tier and each sub-NLZ state, records
the average true cardinality at which buckets occupy that state.  Outputs
per-state correction factors (CFs) in log2 space:

    CF(state) = log2( mean_cardinality_in_state / mean_cardinality_in_tier )

These CFs are the values stored in StateTable.java (CF_HISTORY_1, CF_HISTORY_2,
CF_HISTORY_3, CF_MANTISSA_2, CF_LUCK_*, CF_HISTMANT_*, etc.) and used by
CardStats phase 2a to build tierMult for the Mean/HMean estimators:

    tierMult = 2^( -(cf + CF_OFFSET) ) * invTermCF

Negative CF means the state's average cardinality is BELOW the tier average
(bucket arrived at this tier quickly, so the estimate should be lowered).
Positive CF means the state's average cardinality is ABOVE the tier average
(bucket has been in this tier a long time, so the estimate should be raised).

Workflow:
  1. Run mantissacompare.sh to generate per-state CFs.
  2. Copy the CF arrays from stdout into StateTable.java.
  3. Pick 'tier 8' values for steady-state arrays (CF_HISTORY_N).
  4. Pick tiers 0 through (N-1) for per-tier arrays (CF_HISTORY_N_TIERS).
  5. Run ddlcalibrate.sh with cf=f to generate a v5 CF table with the new CFs.
  6. Run ddlcalibrate.sh with cf=t to verify accuracy.

Usage:  mantissacompare.sh mode=history bits=2

Modes (mode=):
  history       N-bit history shift register (PLL16c/UDLL6 style).
                Tracks whether sub-tiers below the current NLZ were observed.
                bits=1: 2 states.  bits=2: 4 states.  bits=3: 8 states.
  mantissa      N-bit inverted mantissa (fractional NLZ extension).
                bits=2: 4 states.
  andtissa      N-bit AND-tissa (AND of hash with shifted self).
                bits=2: 4 states.
  nlz2          N-bit secondary NLZ (NLZ of remaining bits after primary NLZ).
                bits=2: 4 states.
  luck          N-bit luck gap (difference between best and second-best NLZ).
                bits=1: 2 states.  bits=2: 4 states.  bits=3: 8 states.
  histmant      Combined history + mantissa.  Use hbits= and mbits= instead
                of bits=.  Total states = 2^(hbits+mbits).
                Example: hbits=2 mbits=2 gives 16 combined states.

Parameters:
inner=32768     Elements per trial (max cardinality simulated per bucket).
outer=131072    Number of independent trials.  More = smoother CFs.
                131072 trials at inner=32768 takes ~60 seconds.
maxtier=11      Highest NLZ tier to record.  Tier 8 is the standard reference
                for steady-state CFs.  Tiers above ~11 may have too few
                samples for reliable statistics.
bits=2          State bits for single-mode operation.  Ignored in histmant mode.
hbits=2         History bits for histmant mode.
mbits=2         Mantissa bits for histmant mode.

Output format (stdout):
  Header:  Tier  Total  P(0)  CF(0)  P(1)  CF(1)  ...
  Per-tier rows with observation counts, state probabilities, and CFs.
  Final line: weighted steady-state CFs across tiers 3-maxtier.

  P(s) = fraction of observations in state s at this tier.
  CF(s) = log2(stateAvg / tierAvg) = additive NLZ correction for state s.
  N/A = too few observations (<10) for reliable CF.

Tables generated from this tool:
  StateTable.CF_HISTORY_1       mode=history bits=1, tier 8 CFs
  StateTable.CF_HISTORY_2       mode=history bits=2, tier 8 CFs
  StateTable.CF_HISTORY_3       mode=history bits=3, tier 8 CFs
  StateTable.CF_HISTORY_1_TIERS mode=history bits=1, tiers 0-1
  StateTable.CF_HISTORY_2_TIERS mode=history bits=2, tiers 0-2
  StateTable.CF_HISTORY_3_TIERS mode=history bits=3, tiers 0-3
  StateTable.CF_MANTISSA_2      mode=mantissa bits=2, tier 8 CFs
  StateTable.CF_ANDTISSA_2      mode=andtissa bits=2, tier 8 CFs
  StateTable.CF_NLZ2_2          mode=nlz2 bits=2, tier 8 CFs
  StateTable.CF_LUCK_1          mode=luck bits=1, tier 8 CFs
  StateTable.CF_LUCK_2          mode=luck bits=2, tier 8 CFs
  StateTable.CF_LUCK_3          mode=luck bits=3, tier 8 CFs
  StateTable.CF_LUCK_2_TIERS    mode=luck bits=2, tiers 0-2
  StateTable.CF_HISTMANT_2H2M   mode=histmant hbits=2 mbits=2, tier 8
  StateTable.CF_HISTMANT_1H3M   mode=histmant hbits=1 mbits=3, tier 8
  StateTable.CF_HISTMANT_3H1M   mode=histmant hbits=3 mbits=1, tier 8
  StateTable.CF_HISTMANT_1H3M_TIERS  mode=histmant hbits=1 mbits=3, tiers 0-1

How the CFs are used at runtime (CardStats phase 2a):
  For each filled bucket with absNlz and history pattern h:
    cf = StateTable.historyOffset(absNlz, hbits, h)
    tierMult = 2^(-(cf + CF_OFFSET)) * invTermCF
    corrDif = dif * tierMult
  Where dif = 2^(63-absNlz) is the bucket's raw contribution to the Mean sum.
  The corrected dif values are summed and fed to the Mean formula.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1g" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.MantissaCompare2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
