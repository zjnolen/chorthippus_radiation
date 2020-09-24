#!/bin/bash
#SBATCH -J pop_miss20
#SBATCH --output=slurm_scripts/pop_miss20.out
#SBATCH --clusters=hugemem
#SBATCH --partition=hugemem_std
#SBATCH --cpus-per-task=10
#SBATCH --mem=200000mb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -t 100:00:00

STARTTIME=$(date +%s)

angsd -b ../input_files/lrz_all.bamlist -ref ../input_files/grasshopperRef.fasta -doMajorMinor 1 -GL 1 -doMaf 1 -doCounts 1 -doAbbababa2 1 -sizeFile pop_all.size -useLast 1 -r chr1: -sites ../input_files/neutral_sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 20 -minInd 67 -setMinDepth 168 -SNP_pval 1e-6 -out output/pop_miss20

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
