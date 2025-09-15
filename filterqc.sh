#!/bin/bash

usage(){
echo "
Written by Shijie Yao
Last modified March 22, 2018

Description: Fastq Filter pipeline

Usage:        filterqc.sh in=<file> out=<dir>

Parameters:
in=file         Specify the input fastq or fastq.gz file
out=dir         The output directory
rqcfilterdata=dir Path to RQCFilterData dir
[qc]            do read qc on filtered output

Please contact Shijie Yao at syao @lbl.gov if you encounter any problems.
"
}


pushd . > /dev/null         # save current dir in directory stack
DIR="${BASH_SOURCE[0]}"     # script file, even by source (vs $0)

while [ -h "$DIR" ]; do     # if a link, follow the link to path
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done

cd "$(dirname "$DIR")"      # obtain the abspath of where the script live
DIR="$(pwd)"

popd > /dev/null            # move back to invocation dir

PYDIR="$DIR/pytools"        # abs path to pytools

set=0                       #?

if [ $# -lt 2 ] || [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then    # -z tells if null
    usage
    exit
fi

fastq=""
out=""
doqc=""
rqcfilterdata=""
parse_arg() {
    IFS='=' read -ra toks <<< "$1"    #Convert string to array
    if [ ${#toks[@]} -gt 1 ]; then
        if [ ${toks[0]} == "in" ]; then
            fastq=${toks[1]}
        elif [ ${toks[0]} == "out" ]; then
            out=${toks[1]}
        elif [ ${toks[0]} == "rqcfilterdata" ]; then
            rqcfilterdata=${toks[1]}
        fi
    elif [ ${#toks[@]} -eq 1 ]; then
        if [ ${toks[0]} == "qc" ]; then
            doqc="doqc"
        fi
    fi
}

parse_arg $1
parse_arg $2
parse_arg $3

if [ -z $fastq ] || [ -z $out ]; then
    usage
    exit
fi

if [ ! -e $fastq ]; then
    echo "ERROR - The file not found : $fastq !!"
    exit
fi

CMD="$PYDIR/filter.py -f $fastq -o $out -rdb $rqcfilterdata --prod-type DNA --skip-blast"
if [ ! -z $doqc ]; then
    CMD="$CMD --qc"
fi
echo "$CMD"
eval $CMD
