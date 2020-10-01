#Runs NGSadmix for multiple replicates, preserving all outputs. Can be run with command `bash replicate_NGSadmix_jobs.sh output_dir beagle_file K nreps` Example `bash replicate_NGSadmix_jobs.sh outputs beagle_p1e-2 9 20`
#output_dir: directory for output files
#beagle_file: beagle file name without .beagle.gz extension
#K: K value of analysis
#nreps: Number of replicates to run

OUTPUT_DIR="$1"

beagle_file="$2"

K=$3

if [ ! -d ${OUTPUT_DIR} ]; then

	mkdir ${OUTPUT_DIR}

fi

SCRIPT_DIR="slurm_scripts"

if [ ! -d ${SCRIPT_DIR} ]; then

	mkdir ${SCRIPT_DIR}

fi

for rep in `eval echo {1..$4}`;

	do

echo "#!/bin/bash
#SBATCH -J ${beagle_file}_K${K}_${rep}
#SBATCH --output=${SCRIPT_DIR}/${beagle_file}_K${K}_${rep}.out
#SBATCH --cpus-per-task=4
#SBATCH -t 48:00:00

STARTTIME="'$(date +"%s")'"

module load angsd

NGSadmix -likes ../beagle/outputs/${beagle_file}.beagle.gz -K ${K} -P 4 -o ${OUTPUT_DIR}/${beagle_file}_K${K}_${rep} -minMaf 0.05

"'ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."' > ${SCRIPT_DIR}/${beagle_file}_K${K}_${rep}.sh

sbatch ${SCRIPT_DIR}/${beagle_file}_K${K}_${rep}.sh

done
