#!/bin/bash

usage(){
	echo "Converts BAM format to SAM format."
	echo "Usage: bamlinestreamer.sh input.bam output.sam"
	echo ""
	exit 1
}

if [ $# -lt 2 ]; then
	usage
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

java -ea -Xmx8g -cp "$DIR/current" stream.Bam2Sam "$@"
