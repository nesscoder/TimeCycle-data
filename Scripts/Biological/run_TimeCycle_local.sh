#!/bin/bash
cd ~/Desktop/TimeCycle-data/Scripts/Biological

while IFS=$'\t' read -r fileName fileIn fileOut timepoints reps maxLag; do
	echo "$fileName submitted on `date`"
	Rscript quest_script_main_TimeCycle.R $fileName $fileIn $fileOut $timepoints $reps $maxLag
	sleep 5m
done <"masterFilePath_timeCycle.txt"

echo "DONE"
exit
