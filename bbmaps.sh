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
  search uses maxindel=50/maxindel2=100 by default. Unmapped single reads and
  paired reads with a strong half-read error asymmetry retry at
  retrymaxindel=16000/retrymaxindel2=32000; retryminmapq=0 accepts every mapped
  retry, while a higher value requires that minimum MAPQ before selecting it.
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
Run bbmap.sh -h for the full flag reference.
Java SIMD is detected by the standard BBTools launcher setup.
"
}

addNeuralMapqResources(){
	NEURAL_DEFAULT_ARGS=()
	local enabled=false pairEnabled=false autoEnabled=true controlSeen=false
	local hasNet=false hasLarge=false hasSmall=false
	local hasPairNet=false hasPairLarge=false hasPairSmall=false hasCaps=false
	local arg key value
	for arg in "$@"; do
		key="${arg%%=*}";key="${key,,}"
		if [[ "$arg" == *=* ]];then value="${arg#*=}";else value="";fi
		value="${value,,}"
		case "$key" in
			neuralmapq)
				controlSeen=true;autoEnabled=false
				case "$value" in ''|t|true|1|yes) enabled=true;; *) enabled=false;; esac
				;;
			neuralmapqpair)
				controlSeen=true;autoEnabled=false
				case "$value" in ''|t|true|1|yes) pairEnabled=true;; *) pairEnabled=false;; esac
				;;
			mapqmode)
				controlSeen=true;autoEnabled=false
				case "$value" in
					neural) enabled=true;pairEnabled=false;;
					neuralpaired) enabled=false;pairEnabled=true;;
					legacy) enabled=false;pairEnabled=false;;
					neuralauto|auto) enabled=false;pairEnabled=false;autoEnabled=true;;
				esac
				;;
			neuralmapqnet) hasNet=true;;
			neuralmapqlutlarge) hasLarge=true;;
			neuralmapqlutsmall) hasSmall=true;;
			neuralmapqpairnet) hasPairNet=true;;
			neuralmapqpairlutlarge) hasPairLarge=true;;
			neuralmapqpairlutsmall) hasPairSmall=true;;
			neuralmapqcaps) hasCaps=true;;
		esac
		done
	if ! $controlSeen;then NEURAL_DEFAULT_ARGS+=("mapqmode=neuralauto");fi
	if $enabled || $autoEnabled;then
		if ! $hasNet;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqnet=?bbmaps_mapq_single_v2.bbnet")
		fi
		if ! $hasLarge;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqlutlarge=?neural_mapq/bbmaps_mapq_single_v2_large.tsv")
		fi
		if ! $hasSmall;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqlutsmall=?neural_mapq/bbmaps_mapq_single_v2_small.tsv")
		fi
	fi
	if $pairEnabled || $autoEnabled;then
		if ! $hasPairNet;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqpairnet=?bbmaps_mapq_paired_v2.bbnet")
		fi
		if ! $hasPairLarge;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqpairlutlarge=?neural_mapq/bbmaps_mapq_paired_v2_large.tsv")
		fi
		if ! $hasPairSmall;then
			NEURAL_DEFAULT_ARGS+=("neuralmapqpairlutsmall=?neural_mapq/bbmaps_mapq_paired_v2_small.tsv")
		fi
	fi
	if { $enabled || $pairEnabled || $autoEnabled; } && ! $hasCaps;then
		NEURAL_DEFAULT_ARGS+=("neuralmapqcaps=?neural_mapq/bbmaps_mapq_v2_caps.tsv")
	fi
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMapS build=1 overwrite=true fastareadlen=500 ${NEURAL_DEFAULT_ARGS[*]} $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" align2.BBMapS build=1 overwrite=true fastareadlen=500 "${NEURAL_DEFAULT_ARGS[@]}" "$@"
}

resolveSymlinks
addNeuralMapqResources "$@"
setEnv "$@"
launch "$@"
