#!/bin/bash

usage(){
echo "
BBMapS — BBMap alignment with the Streamer/Writer interface.

Single-ended:  bbmapS.sh ref=reference.fa in=reads.fq out=mapped.sam
Paired-end:    bbmapS.sh ref=reference.fa in=R1.fq in2=R2.fq out=mapped.sam
Index only:    bbmapS.sh ref=reference.fa path=index
Reuse index:   bbmapS.sh in=reads.fq out=mapped.sam path=index
Split reads:   bbsplitS.sh ref_a=a.fa ref_b=b.fa in=reads.fq basename=out_%.fq

Same flag surface as bbmap.sh (build=, in=, in2=, ref=, t=, out=, etc.).
Optional short-indel acceleration: quantumonebase=t (default f).
Optional no-MSA speed mode: quantumonly=t (default f; changes scoring/mapping).
Optional adaptive Quantum/MSA hybrid: quantumhybrid=t (default f; uses Quantum
  for supported short-indel sites and defers compressed MSA when evidence allows).
Optional k-mer pseudoalignment: pseudoalign=t (default f; emits a seed-derived
  polycrystalline CIGAR without full MSA or identity/edit filters; intended for
  coverage/counting).
Optional selective max-indel retry: hybridmaxindel=t (default f). The first
  search uses maxindel=50/maxindel2=100 by default. Unmapped single reads retry;
  paired reads retry when no acceptable pair is found or a mostly clean read
  has clustered errors near one end. The retry uses
  retrymaxindel=16000/retrymaxindel2=32000. Paired retries must pass geometry
  checks; retryminmapq=0 adds no MAPQ restriction, while a higher value also
  requires that minimum MAPQ. Hybrid is conservative and can miss alignments
  found by a full maxindel=16000 search; it remains default-off.
Calibrated neural MAPQ is enabled by default. BBMapS automatically selects the
  single or paired V2 model from the actual input stream. Set mapqmode=legacy
  (or neuralmapq=f/neuralmapqpair=f) to restore legacy MAPQ. Automatic neural
  MAPQ yields to feature-export, match=f, perfect/semiperfect, Quantum, and
  pseudoalignment modes. V2 supports
  primary single-end standard alignments of 50-250 bp. References <=20 Mb use
  the small-reference calibration; references >=1 Gb use the large-reference
  calibration. Length-specific evidence caps range Q26-Q41; unsupported lengths
  and intermediate references retain legacy MAPQ.
  The launcher supplies frozen models and calibration data through BBTools'
  standard ? resource lookup.
  Advanced users may override neuralmapqnet=, neuralmapqlutlarge=, and
  neuralmapqlutsmall=, and neuralmapqcaps= explicitly.
Optional explicit paired selection: neuralmapqpair=t. Paired
  V2 emits one MAPQ per mapped primary mate of 50-250 bp using both ends and
  pair geometry. References <=20 Mb use the small-reference calibration;
  references >=1 Gb use the large-reference calibration. Length-specific caps
  range Q31-Q45; unsupported lengths and intermediate references retain legacy
  MAPQ. The launcher supplies the paired V2 resources.
  Advanced users may override neuralmapqpairnet=,
  neuralmapqpairlutlarge=, and neuralmapqpairlutsmall= explicitly.
BBMapS also computes calibrated strict-endpoint MAPQ for supported primary
  alignments. QL:i stores loose-placement MAPQ and QS:i stores exact reference,
  strand, start, and stop MAPQ. The primary MAPQ field uses QL by default;
  strictmapq=t swaps it to QS while retaining both tags. Unsupported lengths
  and intermediate reference sizes retain legacy MAPQ and omit QL/QS.
Run bbmap.sh -h for the full flag reference.
Java SIMD is detected by the standard BBTools launcher setup.
BBMapS automatically uses horizontal SIMD for limited and unlimited alignment when supported.
  simd=f or an older Java runtime retains the original scalar alignment path.
  The slower experimental diagonal unlimited fill is never selected automatically.
  To try it: JAVA_TOOL_OPTIONS=-Dbbmap3.horizontalUnlimited=diagonal bbmapS.sh ...
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMapS build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" align2.BBMapS build=1 overwrite=true fastareadlen=500 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
