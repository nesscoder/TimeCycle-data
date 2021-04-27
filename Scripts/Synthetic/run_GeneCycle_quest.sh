#!/bin/bash

while IFS=$'\t' read fileName fileIn fileOut interval length reps
do
	here=`pwd`
	JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH -A p30673
#SBATCH --partition=long
#SBATCH --time=120:00:00
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=elanness-cohn2017@u.northwestern.edu
#SBATCH -J $fileName
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --mem=5G

# unload modules that may have been loaded when job was submitted
module purge all

# load the module you want to use 
module load R/3.5.1

# By default all file paths are relative to the directory where you submitted the job.
cd /home/emn6548/TimeCycle/Scripts/Synthetic/

Rscript quest_script_main_GeneCycle.R $fileName $fileIn $fileOut $interval $length $reps
EOJ
`
# print out the job id for reference later
echo "JobID = ${JOB} for parameters $fileName submitted on `date`"
done < masterFilePath_GeneCycle.txt
exit

# make this file executable and then run from the command line
# chmod u+x run_GeneCycle_quest.sh
# ./run_GeneCycle_quest.sh
