NGSadmix
================

-   [NGSadmix in ANGSD](#ngsadmix-in-angsd)
-   [Input files - Making a Beagle file](#input-files---making-a-beagle-file)
    -   [Beagle SLURM script](#beagle-slurm-script)
    -   [Beagle command line options](#beagle-command-line-options)
-   [Running NGSadmix](#running-ngsadmix)
    -   [Replicating and selecting the highest likelihood](#replicating-and-selecting-the-highest-likelihood)
-   [Plotting Results](#plotting-results)
-   [References](#references)

NGSadmix in ANGSD
-----------------

[ANGSD includes NGSadmix](http://www.popgen.dk/angsd/index.php/NGSadmix), a command line software used to estimate population structure from next generation sequence data (Skotte, Korneliussen, and Albrechtsen 2013). This is done by grouping individuals into K clusters, based on Hardy-Weinberg equilibrium. The software is similar to Structure, but can instead use genotype likelihood data as input. This allows for uncertainty to be included in the analysis that takes into account sequencing errors, coverage, and more. This is especially useful for transcriptome data, as coverage varies considerably depending on gene expression. Because both ANGSD and NGSadmix are based on genotype likelihoods, both are packaged together in current releases of ANGSD.

Input files - Making a Beagle file
----------------------------------

NGSadmix uses a Beagle file, which includes genotype likelihood data for a group of individuals. We have uploaded our Beagle files used in this study to the [Dryad](ttps://doi.org/10.5061/dryad.pzgmsbchj) as a tar archive (`beagle_chorthippus_p1e-2.tar.gz`) We note the process for building this file below.

This file was built with ANGSD and requires a bamlist pointing to the bam files. Note that the order of files in the bamlist will correspond to the order of individuals in the final output.

The Beagle file is then built by using ANGSD with the option `-doGlsf 2`. Beagle files take a long time to build and are very large, so it is best to run them as SLURM scripts on the cluster. Below is an example, [`beagle_p1e-2.sh`](beagle_p1e-2.sh):

#### Beagle SLURM script

    #!/bin/bash
    #SBATCH -J beagle_p1e-2
    #SBATCH --output=beagle_p1e-2.out
    #SBATCH --cpus-per-task=6
    #SBATCH --mem=150000mb
    #SBATCH -t 168:00:00

    STARTTIME=$(date +"%s")

    angsd -b ../lrz_all.bamlist -ref grasshopperRef.fasta -doMajorMinor 1 -GL 1 -doGlf 2 -SNP_pval 1e-2 -doMaf 1 -nThreads 2 -r chr1: -sites neutral_sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 20 -doCounts 1 -doPost 2 -doGeno 32 -minInd 67 -setMinDepth 168 -out ../outputs/p1e-2

    ENDTIME=$(date +%s)
    TIMESPEND=$(($ENDTIME - $STARTTIME))
    ((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
    timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
    echo "Took $timestamp hours:minutes:seconds to complete..."

It is important to allot enough memory for the process, here we reserve 150000 MB for a beagle file of all 84 *Chorthippus* and *Pseudochorthippus* individuals. Depending on the cluster you may need to move to a large memory partition or reserve more cores to reserve such a large amount of memory.

#### Beagle command line options

See [General ANGSD Options](../#general-options) for options not mentioned below.

<table>
<colgroup>
<col width="46%" />
<col width="54%" />
</colgroup>
<thead>
<tr class="header">
<th>Option</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>-doGlf 2</code></td>
<td>Dumps genotype log likelihoods to a file, 2 refers to beagle format. <a href="http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods">Source</a></td>
</tr>
<tr class="even">
<td><code>-nThreads 2</code></td>
<td>Sets the number of threads the analysis can run on, this will speed up the analysis, but high values can be unstable.</td>
</tr>
<tr class="odd">
<td><code>-doPost 2</code></td>
<td>Calculates posterior genotype probability, 2 means a uniform prior will be used. <a href="http://www.popgen.dk/angsd/index.php/Genotype_calling">Source</a></td>
</tr>
<tr class="even">
<td><code>-doGeno 32</code></td>
<td>Determines how output will be printed to file, 32 dumps the posterior probabilities in binary. <a href="http://www.popgen.dk/angsd/index.php/Genotype_calling">Source</a></td>
</tr>
<tr class="odd">
<td><code>-out ../outputs/p1e-2</code></td>
<td>Location and base name of output file. When producing multiple files with different options, be sure to include variable options in the file name.</td>
</tr>
</tbody>
</table>

Running NGSadmix
----------------

#### Replicating and selecting the highest likelihood

In order to obtain the highest likelihood clustering of individuals, it is important to run the analysis multiple times. We replicated this analysis 50 times per K cluster, and selected the highest likelihood replicate to report.

For this, we have created a [bash file](replicate_NGSadmix_jobs.sh) that submits a SLURM job with a unique output for each replicate. This can be run from the command line with the following command:

`bash replicate_NGSadmix_jobs.sh output_dir beagle_file K nreps`

where:

| Argument     |                                                   |
|--------------|---------------------------------------------------|
| output\_dir  | Name of the output directory                      |
| beagle\_file | Beagle file name without the .beagle.gz extension |
| K            | Number of clusters in the analysis                |
| nreps        | Total number of replicates to perform             |

For example, to run 50 replicates of a K9 analysis on the p1e-2.beagle.gz file, you would use:

`bash replicate_NGSadmix_jobs.sh outputs p1e-2 9 50`

Plotting Results
----------------

The final results can be plotted with R, however the script needs to be heavily edited based on the bamlist of individuals. For our data, we have a general R script: [plotNGSadmix.R](plotNGSadmix.R) which is used to plot K2-10 for SNP p-value 1e-2 using 50 replicates, which is the results we report.

The main points to edit when changing bamlists are colors, mainlines, dotlines, and spacevec, which build the color palette of the figure and the dividers between populations and individuals. Colors will need to be set per analysis, as the clusters are independently determined (cluster 1 in a K2 analysis may not contain the same group of individuals as cluster 1 in a K3 analysis). In addition, beagle file names will need to be edited, as this script is built around simply changing the p-value in the file name.

References
==========

Skotte, Line, Thorfinn Sand Korneliussen, and Anders Albrechtsen. 2013. “Estimating Individual Admixture Proportions from Next Generation Sequencing Data.” *Genetics* 195 (3): 693–702. doi:[10.1534/genetics.113.154138](https://doi.org/10.1534/genetics.113.154138).
