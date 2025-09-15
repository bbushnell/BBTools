#!/bin/bash

#usage(){
#	echo "CalcMem v1.19"
#	echo "Written by Brian Bushnell, Doug Jacobsen, Alex Copeland, Bryce Foster"
#	echo "Calculates available memory in megabytes"
#	echo "Last modified November 15, 2024"
#}

#Also parses other Java flags
function parseXmx () {
	
	local setxmx=0
	local setxms=0
	
	SIMD=""
	
	#Creating an array of the args forces expansion of environment variables before iteration
	local ALIST=("$@") ;
	for arg in "${ALIST[@]}"
	do
		if [[ "$arg" = "Xmx="* ]] || [[ "$arg" = "xmx="* ]]; then
			z="-Xmx"${arg:4}
			setxmx=1
		elif [[ "$arg" = "-Xmx="* ]] || [[ "$arg" = "-xmx="* ]]; then
			z="-Xmx"${arg:5}
			setxmx=1
		elif [[ "$arg" = "-Xmx"* ]] || [[ "$arg" = "-xmx"* ]]; then
			#z="$arg"
			z="-X"${arg:2}
			setxmx=1
		elif [[ "$arg" = "Xmx"* ]] || [[ "$arg" = "xmx"* ]]; then
			#z="-$arg"
			z="-X"${arg:1}
			setxmx=1
		elif [[ "$arg" = "-Xms"* ]]; then
			z2="$arg"
			setxms=1
		elif [[ "$arg" = "Xms"* ]]; then
			z2="-$arg"
			setxms=1
		elif [[ "$arg" = "-da" ]] || [[ "$arg" = "-ea" ]]; then
			EA="$arg"
		elif [[ "$arg" = "da" ]] || [[ "$arg" = "ea" ]]; then
			EA="-$arg"
		elif [[ "$arg" = "ExitOnOutOfMemoryError" ]] || [[ "$arg" = "exitonoutofmemoryerror" ]] || [[ "$arg" = "eoom" ]]; then
			EOOM="-XX:+ExitOnOutOfMemoryError"
		elif [[ "$arg" = "-ExitOnOutOfMemoryError" ]] || [[ "$arg" = "-exitonoutofmemoryerror" ]] || [[ "$arg" = "-eoom" ]]; then
			EOOM="-XX:+ExitOnOutOfMemoryError"
		elif [[ "$arg" = "json" ]] || [[ "$arg" = "json=t" ]] || [[ "$arg" = "json=true" ]] || [[ "$arg" = "format=json" ]]; then
			json=1
		elif [[ "$arg" = "silent" ]] || [[ "$arg" = "silent=t" ]] || [[ "$arg" = "silent=true" ]]; then
			silent=1
		elif [[ "$arg" = "simd" ]] || [[ "$arg" = "SIMD" ]] || [[ "$arg" = "simd=t" ]] || [[ "$arg" = "simd=true" ]]; then
			SIMD="--add-modules jdk.incubator.vector"
		fi
	done
	
	if [[ "$setxmx" = "1" ]] && [[ "$setxms" = "0" ]]; then
		local substring=`echo $z| cut -d'x' -f 2`
		z2="-Xms$substring"
		setxms=1
	elif [[ "$setxmx" = "0" ]] && [[ "$setxms" = "1" ]]; then
		local substring=`echo $z2| cut -d's' -f 2`
		z="-Xmx$substring"
		setxmx=1
	fi
	
	set=$setxmx
	
}

function setEnvironment(){

	EA="-ea"
	EOOM=""

	if [[ "$SHIFTER_RUNTIME" = "1" ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [ -z "$EC2_HOME" ]; then
		#Let's assume this is the AWS taxonomy server...
		PATH=/test1/binaries/bgzip:$PATH
		PATH=/test1/binaries/lbzip2/bin:$PATH
		PATH=/test1/binaries/sambamba:$PATH
		#PATH=/test1/binaries/java/jdk-11.0.2/bin:$PATH
		PATH=/test1/binaries/pigz2/pigz-2.4:$PATH
	elif [ -z "$NERSC_HOST" ]; then
		#Not NERSC; do nothing
		:
	elif [ -z "$NERSC_HOST" ]; then
		#Not NERSC; do nothing
		:
	else
		PATH=/global/cfs/cdirs/bbtools/bgzip:$PATH
		PATH=/global/cfs/cdirs/bbtools/lbzip2/bin:$PATH
                PATH=/global/cfs/cdirs/bbtools/samtools116/samtools-1.16.1:$PATH
		#PATH=/global/projectb/sandbox/gaag/bbtools/sambamba:$PATH
		PATH=/global/cfs/cdirs/bbtools/java/jdk-17/bin:$PATH
		PATH=/global/cfs/cdirs/bbtools/pigz2/pigz-2.4:$PATH
	fi
}


function freeRam(){
	RAM=0;

	#Memory is in kilobytes.
	local defaultMem=3200000
	if [ $# -gt 0 ]; then
		defaultMem=$1;
		case $defaultMem in
			*g)
			defaultMem=`echo $defaultMem| cut -d'g' -f 1`
			defaultMem=$(( $defaultMem * $(( 1024 * 1024 )) ))
			;;
			*m)
			defaultMem=`echo $defaultMem| cut -d'm' -f 1`
			defaultMem=$(( $defaultMem * 1024 ))
			;;
			*k)
			defaultMem=`echo $defaultMem| cut -d'k' -f 1`
			;;
		esac
	fi

	local mult=84
	if [ $# -gt 1 ]; then
		mult=$2;
	fi
	
	#echo "mult =    $mult" # percent of memory to allocate
	#echo "default = $defaultMem"
	
	local ulimit=$(ulimit -v)
	ulimit="${ulimit:-0}"
	if [ "$ulimit" = "unlimited" ]; then ulimit=0; fi
	local x=$ulimit
	#echo "x = ${x}" # normally ulimit -v
	
	#local HOSTNAME=`hostname`
	local sge_x=0
	local slurm_x=$(( SLURM_MEM_PER_NODE * 1024 ))

	if [[ $RQCMEM -gt 0 ]]; then
		#echo "branch for manual memory"
		x=$(( RQCMEM * 1024 ));
	elif [ -e /proc/meminfo ]; then
		local vfree=$(cat /proc/meminfo | awk -F: 'BEGIN{total=-1;used=-1} /^CommitLimit:/ { total=$2 }; /^Committed_AS:/ { used=$2 } END{ print (total-used) }')
		local pfree=$(cat /proc/meminfo | awk -F: 'BEGIN{free=-1;cached=-1;buffers=-1} /^MemFree:/ { free=$2 }; /^Cached:/ { cached=$2}; /^Buffers:/ { buffers=$2} END{ print (free+cached+buffers) }')
		
		#echo "vfree =   $vfree"
		#echo "pfree =   $pfree"
		#echo "ulimit =  $ulimit"

		local x2=0;

		
		if [ $vfree -gt 0 ] && [ $pfree -gt 0 ]; then
			if [ $vfree -gt $pfree ]; then x2=$pfree; 
			else x2=$vfree; fi
		elif [ $vfree -gt 0 ]; then x2=$vfree;
		elif [ $pfree -gt 0 ]; then x2=$pfree;
		fi

		#echo $sge_x
		#echo $slurm_x
		#echo $x
		#echo $x2

		# set to SGE_HGR_RAMC or SLURM_MEM_PER_NODE value
		if [ $sge_x -gt 0 ]; then 
			if [ $x2 -gt $sge_x ] || [ $x2 -eq 0 ]; then 
				x=$sge_x;
				x2=$x; 
			fi
		elif [ $slurm_x -gt 0 ]; then
			if [ $x2 -gt $slurm_x ] || [ $x2 -eq 0 ]; then 
				x=$slurm_x;
				x2=$x; 
			fi
		fi
		
		#echo "x = ${x}"
		#echo "x2 = ${x2}"
		#echo $vfree
		#echo $pfree
		
		if [ "$x" = "unlimited" ] || (("$x" > $x2)); then x=$x2; fi
		if [ $x -lt 1 ]; then x=$x2; fi
	fi

#TODO:
#"NERSC_HOST", "perlmutter"

	if [ $x -lt 1 ] || [[ "$HOSTNAME" = "genepool"* ]]; then
		#echo "branch for unknown memory"
		#echo $x
		#echo "ram is unlimited"
		RAM=$((defaultMem/1024))
		echo "Max memory cannot be determined.  Attempting to use $RAM MB." 1>&2
		echo "If this fails, please add the -Xmx flag (e.g. -Xmx24g) to your command, " 1>&2
		echo "or run this program qsubbed or from a qlogin session on Genepool, or set ulimit to an appropriate value." 1>&2
	else
		#echo "branch for known memory"
		#echo "x = ${x}"
		#echo "m = ${mult}"
		
		# available (ram - 500k) * 85% / 1024kb = megs of ram to use
		# not sure where this formula came from
		RAM=$(( ((x-500000)*mult/100)/1024 ))
		#echo $RAM
	fi
	#local z="-Xmx${RAM}m"
	#echo $RAM
	return 0
}

#freeRam "$@"
