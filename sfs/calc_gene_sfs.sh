#!/bin/bash
#SBATCH -J calculate_per_gene_sfs_$taxa
#SBATCH --output=calc_per_gene_sfs_$taxa.out
#SBATCH --cpus-per-task=1
#SBATCH --clusters=mpp2
#SBATCH --time=1-00:00:00

#Builds SFS for each gene within specified species pairs. Requires saf files for each population going in.
#These take a very long time to run. Approx. 24 hours per pair. It might be best to split this up into multiple runs so it doesn't take several days. Replace instances of $taxa with population pairs to be run on. This script is written for the annotation pop1_pop2 when inputting taxa and will automatically split this combination.

#Load list of genes from reference. This is a file that has the positions of all genes on artificial chromosome 1, which we limited our analyses to.
genes=`cat grasshopperRef_chr1genes`

	#Split taxa pair into taxa 1 and 2.
	count=${taxa}
	half=$((count-1))
	half=$((half/2))

	first=$(expr substr ${taxa} 1 $half)
	second=$(expr substr ${taxa} $((half+2)) $count)

	echo "Building per gene 2D SFS from the following taxa: ${first} and ${second}"

	if [ ! -d sfs/per_gene_sfs/${taxa} ]; then

		mkdir -p sfs/per_gene_sfs/${taxa}

	fi

	if [ "$second" == "crub_a" ]; then

		second="crub"

	fi

	#Loop over all genes
	for gene in $genes; do

		gene_name=`echo "$gene" | cut -c 6-`

		if [ -f sfs/per_gene_sfs/${taxa}/2dsfs_${taxa}_p1_fold0_${gene_name}.sfs ]; then

			echo "Skipping ${gene}, already analyzed..."

		else

			echo ""
			echo "Building SFS for region: $gene"
			echo ""

			echo $gene_name

			if [ $taxa == 'cppar_cpery' ]; then

				realSFS saf/${first}_p1_fold0_ancmol.saf.idx saf/${second}_p1_fold0_ancmol.saf.idx -r $gene -P 1 > sfs/per_gene_sfs/${taxa}/2dsfs_${taxa}_p1_fold0_${gene_name}.sfs

			else

				realSFS saf/${first}_p1_fold0.saf.idx saf/${second}_p1_fold0.saf.idx -r $gene -P 1 > sfs/per_gene_sfs/${taxa}/2dsfs_${taxa}_p1_fold0_${gene_name}.sfs

			fi

		fi


	done


	#Delete any SFS without data
	find sfs/per_gene_sfs/$taxa -size 0 -delete

	echo "All gene SFS complete and empty SFS deleted! Directory is ready for making summed SFS and bootstrap SFS"

#done
