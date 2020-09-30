#!/bin/bash
#SBATCH -J calc_glob_fst
#SBATCH --output=calc_glob_fst.out
#SBATCH --partition=mpp2_batch
#SBATCH --clusters=mpp2
#SBATCH --cpus-per-task=1
#SBATCH -t 2:00:00

for pair in "cbig_b_cbru_b" "cbig_b_cbru_s" "cbig_e_cbru_b" "cbig_e_cbru_s" "cmol_b_cbru_b" "cmol_b_cbru_s" "cmol_e_cbru_b" "cmol_e_cbru_s" "cmol_e_cbig_e" "cmol_b_cbig_b" "cmol_e_cbig_b" "cmol_b_cbig_e"

	do

		first=$(expr substr ${pair} 1 6)
		second=$(expr substr ${pair} 8 13)

		echo "Calculating global Fst from 2D SFS for the following taxa: ${first} and ${second}"


		realSFS fst index saf/${first}_p1_fold0.saf.idx saf/${second}_p1_fold0.saf.idx -sfs sfs/2dsfs_${first}_${second}_p1_fold0.ml -fstout fst/glob_${first}_${second}_p1_fold0

done
