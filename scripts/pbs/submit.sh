#!/bin/bash

JOBSCRIPT=$1
GENOME=$2

for exp in ???? 
do 
	pushd $exp
	for sample in ????
	do 
		pushd $sample
		for replicate in ?
		do 
			pushd $replicate
			qsub $JOBSCRIPT -N ${exp}/${sample}/${replicate} -F $GENOME
			sleep 2
			popd
			
		done
		popd
	done
	popd
done
