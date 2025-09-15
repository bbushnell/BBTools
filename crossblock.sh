#!/bin/bash

#For usage information, please see decontaminate.sh

function crossblock(){
	CMD="decontaminate.sh $@"
	eval $CMD
}

crossblock "$@"

