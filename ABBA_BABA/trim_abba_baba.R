#This function is used to trim the output of angsd's doABBABABA test down to only the topologies following the known phylogenetic tree and save them as a csv.
#angsd calculates a D-statistic for every possible combination of populations input.
#Not all combinations are valid, as for the test to be informative, the topology tested must match the known phylogenetic tree
#This script takes the tables produced by angsd as an input, along with the species tree described below to remove any outputs that violate the species tree

#IMPORTANT: Be sure to change the species tree if needed before running the script. As of now, it is set to the phylogeny for Chorthippus spp. produced by Burçin Yildirim

library(data.table)
library(ape)

args <- commandArgs(trailingOnly=TRUE)

#Known phylogenetic tree for populations
species_tree <- read.tree(text = "(((((Cbru_S,Cbru_B),Crub),(Cbig_B,Cbig_E)),(Cmol_E,Cmol_B)),Cpar);")

message("This is your selected species tree.")
species_tree

dataset <- fread(args[1],header = TRUE)
filename <- tools::file_path_sans_ext(args[1])

topologies <- dataset[,9:12]

checkassumption <- function(x) {

  tips <- unlist(x)
  taxa <- species_tree$tip.label

  test.tree <- read.tree(text = paste("(((",tips[1],",",tips[2],"),",tips[3],"),",tips[4],");",sep =""))
  true.tree <- drop.tip(species_tree,taxa[!taxa %in% tips])

  all.equal.phylo(test.tree,true.tree)

}

truetopologies <- apply(topologies, 1, checkassumption)

write.table(dataset[truetopologies,], paste(filename,".trimmed.csv",sep = ""), append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
