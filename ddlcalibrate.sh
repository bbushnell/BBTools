#!/bin/bash

usage(){
echo "
Written by Brian Bushnell, Chloe, and Eru
Last modified April 6, 2026

Description:  Calibrates DynamicDemiLog cardinality estimators by feeding random
longs to DDL instances with varied seeds, tracking true cardinality via PRNG counter,
and reporting per-estimator accuracy statistics at logarithmically spaced cardinality
checkpoints.  Cache-friendly: each thread processes one DDL at a time (create, feed
all values, record stats at thresholds, discard), keeping the DDL's working set in
L1/L2 cache throughout.  Results are merged after all threads complete.

Usage:  ddlcalibrate.sh loglogtype=ddl ddls=1000 buckets=2048

Estimator types (loglogtype= or type=):
  ddl (ddl10)   DynamicDemiLog, 10-bit mantissa.
  ddl2          DynamicDemiLog2, 2-bit mantissa.
  ddl8          DynamicDemiLog8, 8-bit mantissa.
  dll2          DynamicLogLog2, 2-bit registers.
  dll3          DynamicLogLog3, 3-bit registers (1-bit mantissa).
  dll3v2        DynamicLogLog3v2, variant with social promotion.
  dll4          DynamicLogLog4, 4-bit registers (no mantissa).
  dll4m         DynamicLogLog4m, mantissa variant.
  bdll3         BankedDynamicLogLog3, banked 3-bit registers.
  ll6           LogLog6, 6-bit registers (no tier promotion).
  udll6         UltraDynamicLogLog6, 6-bit with 2-bit history.
  pll16c        ProtoLogLog16c, 16-bit with configurable mode bits.
  htb           HyperTwoBits, 2-bit threshold estimator.
  htc           HLLTailCut, tail-cutoff HLL variant.
  ull8          UltraLogLog8, 8-bit with history.
  ertl          ErtlULL, Ertl's improved estimator.

Parameters:
ddls=1000       Number of DDL instances with varied seeds (distributed across threads).
buckets=2048    Buckets per DDL instance. Must be a multiple of 256.
maxmult=10      Stop when trueCard reaches buckets*maxmult.
dupfactor=0     Number of duplicate adds per unique value (0=none).
reportfrac=0.01 Report every this fraction of current cardinality (cascading:
                gives one row per add at low cardinality, one per 10 at 1000, etc).
seed=12345      Master seed for DDL seed generation.
threads=1       Number of parallel simulation threads.
out3=           Output file for v5 CF table (File 3).  Not written if omitted.
out4=           Output file for per-DLC-tier data.  Not written if omitted.
cf=f            Set cf=t to apply existing correction factors during calibration.
                Requires matching *CorrectionFactor.tsv in resources/.
formulas=f      Use closed-form CF formulas instead of tables (for whitelisted types).
clamp=t         Clamp estimates to added count (clamptoadded).

Overflow control (DLL2/DLL3/BDLL3 only):
io=t            Ignore overflow buckets in estimation (ignoreoverflow).
co=t            Apply overflow correction (correctoverflow).
ep=t            Early promotion (advance tier when all buckets nonzero).

DLC tuning:
dlcalpha=0.25   DLC logspace alpha.
dlcblendlo=     DLC blend low threshold (fraction of buckets).
dlcblendhi=     DLC blend high threshold (fraction of buckets).

Formula control:
meancfformula=f     Use Mean CF formula (per-class).
hccfformula=f       Use HC CF formula (UDLL6 only).
sbsformula=f        Use SBS formula instead of SBS table.
useformulas=f       Enable all available formulas (alias: formulas).

ProtoLogLog16 mode bits (PLL16c only):
hbits=2         History bits (1, 2, or 3).
lbits=0         Luck bits.
mbits=0         Mantissa bits.
pllmode=        Mode: mantissa, andtissa, nlz2, history, luck, histmant, none.

File 1 output columns (written after all threads complete):
TrueCard        Ground-truth distinct count.
Occupancy       Fraction of DDL buckets with any data (averaged across all DDLs).
*_err           Signed relative error: (estimate - true) / true.  Zero = unbiased.
*_abs           Mean absolute relative error.  Zero = perfectly accurate.
*_std           Population stdev of relative error across DDL instances.  Zero = no variance.

File 3 output columns (out3=, v4 CF table):
TrueCard        Integer ground-truth cardinality (key for CF lookup).
*_cf            Correction factor per estimator: 1/(1+avgErr).  Apply at runtime.
                9 columns: Mean, HMean, HMeanM, GMean, Hybrid, DLC, DLCBest, DLC3B, DThHyb.
                Copy to resources/*CorrectionFactor.tsv for runtime use with cf=t.

Estimators reported:
Mean            Arithmetic mean of bucket values, with correction.
HMean           HLL-style harmonic mean, with correction.
HMeanM          Harmonic mean with mantissa bits (DDL only), with correction.
GMean           Geometric mean proxy, with correction.
HLL             HyperLogLog estimator.
LC              LinearCounting estimator (accurate only at low occupancy).
Hybrid          Blend of LC and Mean/HMeanM based on occupancy.
HybDLC50        Hybrid using DLC for zone detection (50% blend threshold).
DThHyb          Type-aware DLC-threshold hybrid (LC->Mean->HMeanM blend).
LCmin           Tier-compensated LC: lcPure * 2^minZeros.
RawDup          Raw duplicate fraction (diagnostic only).
DLC             DynamicLogLog-Corrected estimate (logspace 0.25).
DLC3B           DLC with 3-bucket smoothing.
DLCBest         Best single-tier DLC estimate.
HybDLC          Hybrid using DLC for zone detection.
DLC0..DLCn      Per-tier raw DLC estimates (diagnostic).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.DDLCalibrationDriver2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
