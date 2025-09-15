#!/bin/bash
#calcmem

#function usage(){
#	echo "CalcMem v1.01"
#	echo "Written by Brian Bushnell, Doug Jacobsen, Alex Copeland"
#	echo "Calculates available memory in megabytes"
#	echo "Last modified April 25, 2014"
#}

function parseXmx () {
	local setxmx=0
	local setxms=0
	
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
			setxmx=1
		elif [[ "$arg" == Xmx* ]]; then
			z="-$arg"
			setxmx=1
		elif [[ "$arg" == -Xms* ]]; then
			z2="$arg"
			setxms=1
		elif [[ "$arg" == Xms* ]]; then
			z2="-$arg"
			setxms=1
		elif [[ "$arg" == -da ]] || [[ "$arg" == -ea ]]; then
			EA="$arg"
		fi
	done
	
	if [[ $setxmx == 1 ]] && [[ $setxms == 0 ]]; then
		local substring=`echo $z| cut -d'x' -f 2`
		z2="-Xms$substring"
		setxms=1
	elif [[ $setxmx == 0 ]] && [[ $setxms == 1 ]]; then
		local substring=`echo $z2| cut -d's' -f 2`
		z="-Xmx$substring"
		setxmx=1
	fi
	
	set=$setxmx
}


RAM=0;

function freeRam(){
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
	
	#echo "mult =    $mult"
	#echo "default = $defaultMem"
	
	local ulimit=$(ulimit -v)
	local x=$ulimit
	
	if [ -e /proc/meminfo ]; then
		local vfree=$(cat /proc/meminfo | awk -F: 'BEGIN{total=-1;used=-1} /^CommitLimit:/ { total=$2 }; /^Committed_AS:/ { used=$2 } END{ print (total-used) }')
		local pfree=$(cat /proc/meminfo | awk -F: 'BEGIN{free=-1;cached=-1;buffers=-1} /^MemFree:/ { free=$2 }; /^Cached:/ { cached=$2}; /^Buffers:/ { buffers=$2} END{ print (free+cached+buffers) }')
		
		echo "vfree =   $vfree"
		echo "pfree =   $pfree"
		echo "ulimit =  $ulimit"

		local x2=0;
		
		if [ $vfree -gt 0 ] && [ $pfree -gt 0 ]; then
			if [ $vfree -gt $pfree ]; then x2=$pfree; 
			else x2=$vfree; fi
		elif [ $vfree -gt 0 ]; then x2=$vfree;
		elif [ $pfree -gt 0 ]; then x2=$pfree;
		fi
		
		echo $x
		echo $x2
		echo $vfree
		echo $pfree
		
		if [ "$x" = "unlimited" ] || (($x > $x2)); then x=$x2; fi
		
	fi
	
	#echo "x=$x"
	local HOSTNAME=`hostname`
	if [ $x -lt 1 ] || [[ $HOSTNAME == genepool* ]]; then
		#echo "hello 2"
		#echo $x
		#echo "ram is unlimited"
		RAM=$((defaultMem/1024))
		echo "Max memory cannot be determined.  Attempting to use $RAM MB." 1>&2
		echo "If this fails, please set ulimit or run this program qsubbed or from a qlogin session on Genepool." 1>&2
	else
		#echo "hello 1"
		#echo $x
		
		#if [ $x -ge 1000000000 ]; then
		#	echo "ram is 1000g+"
		#elif [ $x -ge 500000000 ]; then
		#	echo "ram is 500g+"
		#elif [ $x -ge 250000000 ]; then
		#	echo "ram is 250g+"
		#elif [ $x -ge 144000000 ]; then
		#	echo "ram is 144g+"
		#elif [ $x -ge 120000000 ]; then
		#	echo "ram is 120g+"
		#elif [ $x -ge 40000000 ]; then
		#	echo "ram is 40g+"
		#else
		#	echo "ram is under 40g"
		#fi
		#echo $x
		RAM=$(( ((x-500000)*mult/100)/1024 ))
		#echo $RAM
	fi
	
	echo $set
	echo $z
	echo $z2
	
	#local z="-Xmx${RAM}m"
	return 0
}

parseXmx "$@"
freeRam "$@"
