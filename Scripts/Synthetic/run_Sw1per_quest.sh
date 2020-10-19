#!/bin/bash
while IFS=$'\t' read fileName fileIn fileOut cycles moveAvg
do
	JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH -A p30673
#SBATCH --partition=normal
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=elanness-cohn2017@u.northwestern.edu
#SBATCH -J $fileName
#SBATCH --nodes=1
#SBATCH -n 4
#SBATCH --mem=40G

## set your working directory

## copy necessary files to run folders
# cp ../SW1PerS_v1.m ./
# cp ../run_SW1PerS.m ./
# cp ${fileIn} ${fileOut}

# unload modules that may have been loaded when job was submitted
module purge all

## job commands; run_Sw1pers is the MATLAB .m file, specified without the .m extension

module load matlab/r2018a
matlab -nosplash -nodesktop -singleCompThread -r "data_file_path = '${fileIn}'; out_dir_path = '${fileOut}'; num_cycles = ${cycles}; movingaverage_window = ${moveAvg}; run_SW1PerS"

EOJ`

echo "JobID = ${JOB} for parameters ${filename} submitted on `date`"
done < masterFilePath_sw1per_failed.txt

exit

# make this file executable and then run from the command line
# chmod u+x run_Sw1per_quest.sh
# ./run_Sw1per_quest.sh
