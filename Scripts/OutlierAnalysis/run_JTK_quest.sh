#!/bin/bash

while IFS=$'\t' read fileName fileIn fileOut timepoints reps start end interval
do
    here=`pwd`
    JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH -A p30673
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=elanness-cohn2017@u.northwestern.edu
#SBATCH -J $fileName
#SBATCH --nodes=1
#SBATCH --mem=3G

# unload modules that may have been loaded when job was submitted
module purge all

# load the module you want to use 
module load R/3.5.1

## set your working directory 
cd /home/emn6548/TimeCycle/Scripts/Synthetic/

Rscript quest_script_main_JTK_CYCLE.R $fileName $fileIn $fileOut $timepoints $reps $start $end $interval
EOJ
`
# print out the job id for reference later
echo "JobID = ${JOB} for parameters $fileName submitted on `date`"
done < masterFilePath_JTK_failed.txt
exit

# make this file executable and then run from the command line
# chmod u+x run_JTK_quest.sh
# ./run_JTK_quest.sh
