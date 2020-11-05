Documentation for SFS calculation available [here](https://zjnolen.github.io/chorthippus_radiation/#/sfs).

File description:

- sfs/ - Data SFS for each population comparison.
- sfs/boot_sfs - Bootstrap SFS calculated for each population comparison.
- sfs/per_gene_sfs - SFS for each gene for each population comparison.
- build_gl_boots.R - Used to combine gene SFS into a single whole genome SFS and 100 resambled bootstrap SFS.
- calc_gene_sfs.sh - Used to calculate SFS per gene from SAF files.
- grasshopperRef_chr1genes - List of gene names and positions for single copy genes in the *Pseudochorthippus* reference transcriptome.
