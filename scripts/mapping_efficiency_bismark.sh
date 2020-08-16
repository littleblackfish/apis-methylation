#for f in ????/????/?/results/bismark_summary/bismark_summary_report.txt; do awk '{if (NR>1) { print $3/$2 }}' $f; done|sort|less


#!/bin/bash

JOBSCRIPT=$1

for exp in ???? 
do 
	for sample in $exp/????
	do 
		for replicate in $sample/?
		do 
		
		echo $exp $(basename $sample) $(basename $replicate) $(awk '{if (NR>1) { s+= $3/$2 }} END{print s/(NR-1)}' $replicate/results/bismark_summary/bismark_summary_report.txt )
		done
	done
done
