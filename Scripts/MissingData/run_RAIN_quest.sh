#!/bin/bash

while IFS=$'\t' read fileName fileIn fileOut interval period peak1 peak2 missing
do
    here=`pwd`
    JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH -A p30673
#SBATCH --partition=short
#SBATCH --time=00:15:00
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
cd /projects/p30673/TimeCycleV3/Scripts/MissingData/

Rscript quest_script_main_RAIN.R $fileName $fileIn $fileOut $interval $period $peak1 $peak2 $missing
EOJ
`
# print out the job id for reference later
echo "JobID = ${JOB} for parameters $fileName submitted on `date`"
done < masterFilePath_RAIN.txt
exit

# make this file executable and then run from the command line
# chmod u+x run_RAIN_quest.sh
# ./run_RAIN_quest.sh