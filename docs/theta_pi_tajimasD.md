Calculating per gene Theta, Pi, and Tajima's D
=========================================

## SAF Files

To calculate per gene theta and Tajima's D, first we built a SAF file for each population as described in the [SFS section](sfs.md#creating-saf-files).

## Creating SFS

Next, we build an SFS for each population with the command:

```bash
realSFS $population.saf.idx -P $num_threads > $population.sfs
```

## Calculating Theta and Pi

### Per site

Per site theta is calculated with a single command:

```bash
angsd \
-bam [population].bamlist -ref grasshopperRef.fasta -anc grasshopperRef.fasta \
-r chr1: -sites neutral_sites -doThetas 1 -doSaf 1 \
-pest $population.sfs -GL 1 -baq 1 -remove_bads 1 -uniqueOnly 1 \
-C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 20 \
-setMinDepth $integer -doCounts 1 -minind $integer -out $population
```

where the bamlist is the same as used to produce the population SAF file, and the SFS is the one produced above. `-setMinDepth` and `-minind` define the minimum depth across all individuals for a site to be included and the minimum number of individuals a site must be present in to be included respectively. We set `-setMinDepth` to two times the number of individuals in a population sample and `-minind` to 0.8 times the number of individuals in the population sample.

This provides a file, `$population.thetas.gz`. We then used

```bash
thetaStat make_bed $population.thetas.gz
```

to produce `$population.thetas.idx`. `thetaStat` is a command of ANGSD.

### Per gene

We use a custom Perl script ([`loopTheta.pl`](https://github.com/zjnolen/chorthippus_radiation/blob/master/theta_pi_tajimasD/loopTheta.pl)) that is specific to our system to then calculate per gene theta:

```bash
perl loopTheta.pl grasshopperRef.positions $population.thetas.idx > \
PG_Thetas_${population}_${sample_size}
```

This produces a table, with the genes as rows and the following columns:

- Column 1 and 2 showed the ID and the location of the genes
- Column 3 was the total length of the gene
- Column 4 was the effective length (effective number of sites)
- Column 5 and 6 were the sum of Watterson Theta and pairwise Pi, respectively.

To calculate the Theta and Pi values for each gene, we divided the sum of Watterson Theta and Pi by the effective length and added them to the file as separate columns.

## Calculating Tajima's D

### Per site

Tajima's D can be calculated from the pairwise Pi values calculated with Theta:

```bash
thetaStat do_stat $population.thetas.idx -nChr 1
```

The output of this command is a file with `.pestPG` extension. The `.pestP`G file has 14 columns with the 5 different Theta estimations and 5 neutrality test statistics, Tajima's D being in the 9th column.

### Per gene

As for Theta and Pi we used a custom Perl script ([loopTajimaD.pl](https://github.com/zjnolen/chorthippus_radiation/blob/master/theta_pi_tajimasD/loopTajimaD.pl)) to calculate per gene Tajima's D:

```bash
perl loopTajimaD.pl grasshopperRef.positions $population.thetas.idx $chr_num > \
PG_Tajima_${population}_${sample_size}
```

where `$chr_num` is the chromosome number. Here, it is the number of individuals in the sample times the ploidy (x2 for diploid).

The output has the same format as the per site file, but with a row for each gene.

We plot Theta and Tajima's D as per gene density plots in Figure S11 and 3c respectively. Below is the density plots for each:

![](_images/theta_tajima_dens.svg)

-------------

<div style="text-align: right">Page last updated on {docsify-updated}</div>
