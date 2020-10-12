#install.packages("PopGenome")
library(PopGenome)

# Alignmets that I prepare for IQTree
Alignments_md10<- readData("alignments_md10")

get.individuals(Alignments_md10)
## Set populations pop1: Cpe, pop2: Cpp, pop3: Cbig, pop4: Cmol, pop5: Cbru, pop6: Crub

Alignments_md10_popdefined <- set.populations(Alignments_md10, list(c("Cpe_19esc","Cpe_21esc","Cpe_347bie","Cpe_348bie","Cpe_349bie"),c("Cpp_25gab","Cpp_27gab", "Cpp_7aru", "Cpp_8aru", "Cpp_9aru"),
                                                         c("Cbig_Bf1", "Cbig_Bf2", "Cbig_Bf3", "Cbig_Bm1", "Cbig_Bm2", "Cbig_Bm3", "Cbig_Bm4", "Cbig_Bm5","Cbig_Ef1","Cbig_Ef2","Cbig_Ef3","Cbig_Ef4","Cbig_Ef5", "Cbig_Em1","Cbig_Em2","Cbig_Em3","Cbig_Em4","Cbig_Em5"),
                                                         c("Cmol_Bf1", "Cmol_Bf2", "Cmol_Bf3", "Cmol_Bf4", "Cmol_Bf5" ,"Cmol_Bm1", "Cmol_Bm2", "Cmol_Bm3", "Cmol_Bm4", "Cmol_Bm5","Cmol_Ef1","Cmol_Ef2","Cmol_Ef3","Cmol_Ef4","Cmol_Ef5", "Cmol_Em1","Cmol_Em2","Cmol_Em3","Cmol_Em4","Cmol_Em5"),
                                                         c("Cbru_Bf1", "Cbru_Bf2", "Cbru_Bf3", "Cbru_Bf4", "Cbru_Bf5" ,"Cbru_Bm1", "Cbru_Bm2", "Cbru_Bm3", "Cbru_Bm4", "Cbru_Bm5","Cbru_Sf1", "Cbru_Sf2", "Cbru_Sf3", "Cbru_Sf4", "Cbru_Sf5" ,"Cbru_Sm1", "Cbru_Sm2", "Cbru_Sm3", "Cbru_Sm4", "Cbru_Sm5"),
                                                         c("Crubratibialis_1", "Crubratibialis_2","Crubratibialis_3","Crubratibialis_4","Crubratibialis_5","Crubratibialis_6","Crubratibialis_7","Crubratibialis_8","Crubratibialis_9","Crubratibialis_10","Crubratibialis_12","Crubratibialis_14","Crubratibialis_15","Crubratibialis_16","Crubratibialis_19", "Crubratibialis_20")))

## calculate dxy
# Dxy calculation : haplotype and nucleotide mode :TRUE, prior definition of populations
trial_1<- diversity.stats.between(Alignments_md10_popdefined,subsites = FALSE, haplotype.mode = TRUE, nucleotide.mode = TRUE)
trial_1@populations # Defined populations beforehand: 1: Cpe, 2: Cpp, 3: Cbig, 4: Cmol, 5: Cbru, 6: Crub

trial_1@poppairs
trial_1@nuc.diversity.between
#trial_1@hap.diversity.between

nucdiv_between_1 <- as.data.frame(trial_1@nuc.diversity.between/trial_1@n.sites)
nucdiv_between_1
#hapdiv_between_1 <- as.data.frame(trial_1@hap.diversity.between)

colnames(nucdiv_between_1) <- c("Cpe_Cpp", "Cpe_Cbig","Cpe_Cmol", "Cpe_Cbru", "Cpe_Crub","Cpp_Cbig","Cpp_Cmol", "Cpp_Cbru","Cpp_Crub", "Cbig_Cmol", "Cbig_Cbru", "Cbig_Crub", "Cmol_Cbru", "Cmol_Crub","Cbru_Crub")
#colnames(hapdiv_between_1) <- c("Cpe_Cpp", "Cpe_Cbig","Cpe_Cmol", "Cpe_Cbru", "Cpe_Crub","Cpp_Cbig","Cpp_Cmol", "Cpp_Cbru","Cpp_Crub", "Cbig_Cmol", "Cbig_Cbru", "Cbig_Crub", "Cmol_Cbru", "Cmol_Crub","Cbru_Crub")

#write.table(nucdiv_between_1, file = "nucdiv_between", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
#write.table(hapdiv_between_1, file = "hapdiv_between", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


#library(readr)
# nucdiv_between <- read_delim("nucdiv_between",
#                              "\t", escape_double = FALSE, col_names = TRUE,
#                              trim_ws = TRUE)
# hapdiv_between <- read_delim("hapdiv_between",
#                              "\t", escape_double = FALSE, col_names = TRUE,
#                              trim_ws = TRUE)


genome_average<- colMeans(nucdiv_between_1)
genome_average
names <- c("Cpe_Cpp", "Cpe_Cbig","Cpe_Cmol", "Cpe_Cbru", "Cpe_Crub","Cpp_Cbig","Cpp_Cmol", "Cpp_Cbru","Cpp_Crub", "Cbig_Cmol", "Cbig_Cbru", "Cbig_Crub", "Cmol_Cbru", "Cmol_Crub","Cbru_Crub")
Dxy_genome_average <- data.frame(Names=names, Dxy_genome_average=genome_average)
p<- ggplot(Dxy_genome_average, aes(x=Names, y=Dxy_genome_average))+geom_bar(stat = "identity",width = 0.6)
p+ylab("Dxy genome average")+theme_minimal()+theme(axis.text.x = element_text(angle = 90))

# colMeans(hapdiv_between_1)

## make a distribution plots.  just manipulate the data frame if you want to use ggplot
# pop1: Cpe, pop2: Cpp, pop3: Cbig, pop4: Cmol, pop5: Cbru, pop6: Crub
library(dplyr)

Cpe_Cpp<- nucdiv_between_1[,1, drop=FALSE]
colnames(Cpe_Cpp) <- "Dxy_nucleotide"
Cpe_Cpp$Name <- "Cpe_Cpp"

Cpe_Cbig<- nucdiv_between_1[,2, drop=FALSE]
colnames(Cpe_Cbig) <- "Dxy_nucleotide"
Cpe_Cbig$Name <- "Cpe_Cbig"

Cpe_Cmol<- nucdiv_between_1[,3, drop=FALSE]
colnames(Cpe_Cmol) <- "Dxy_nucleotide"
Cpe_Cmol$Name <- "Cpe_Cmol"

Cpe_Cbru<- nucdiv_between_1[,4, drop=FALSE]
colnames(Cpe_Cbru) <- "Dxy_nucleotide"
Cpe_Cbru$Name <- "Cpe_Cbru"

Cpe_Crub<- nucdiv_between_1[,5, drop=FALSE]
colnames(Cpe_Crub) <- "Dxy_nucleotide"
Cpe_Crub$Name <- "Cpe_Crub"

Cpp_Cbig<- nucdiv_between_1[,6, drop=FALSE]
colnames(Cpp_Cbig) <- "Dxy_nucleotide"
Cpp_Cbig$Name <- "Cpp_Cbig"

Cpp_Cmol<- nucdiv_between_1[,7, drop=FALSE]
colnames(Cpp_Cmol) <- "Dxy_nucleotide"
Cpp_Cmol$Name <- "Cpp_Cmol"

Cpp_Cbru<- nucdiv_between_1[,8, drop=FALSE]
colnames(Cpp_Cbru) <- "Dxy_nucleotide"
Cpp_Cbru$Name <- "Cpp_Cbru"

Cpp_Crub<- nucdiv_between_1[,9, drop=FALSE]
colnames(Cpp_Crub) <- "Dxy_nucleotide"
Cpp_Crub$Name <- "Cpp_Crub"

Cbig_Cmol<- nucdiv_between_1[,10, drop=FALSE]
colnames(Cbig_Cmol) <- "Dxy_nucleotide"
Cbig_Cmol$Name <- "Cbig_Cmol"

Cbig_Cbru<- nucdiv_between_1[,11, drop=FALSE]
colnames(Cbig_Cbru) <- "Dxy_nucleotide"
Cbig_Cbru$Name <- "Cbig_Cbru"

Cbig_Crub<- nucdiv_between_1[,12, drop=FALSE]
colnames(Cbig_Crub) <- "Dxy_nucleotide"
Cbig_Crub$Name <- "Cbig_Crub"

Cmol_Cbru<- nucdiv_between_1[,13, drop=FALSE]
colnames(Cmol_Cbru) <- "Dxy_nucleotide"
Cmol_Cbru$Name <- "Cmol_Cbru"

Cmol_Crub<- nucdiv_between_1[,14, drop=FALSE]
colnames(Cmol_Crub) <- "Dxy_nucleotide"
Cmol_Crub$Name <- "Cmol_Crub"

Cbru_Crub<- nucdiv_between_1[,15, drop=FALSE]
colnames(Cbru_Crub) <- "Dxy_nucleotide"
Cbru_Crub$Name <- "Cbru_Crub"

genes<- as.data.frame(row.names(nucdiv_between_1))
Dxy_nucleotide_all_1 <- rbind(Cpe_Cpp, Cpe_Cbig,Cpe_Cmol, Cpe_Cbru, Cpe_Crub,Cpp_Cbig,Cpp_Cmol, Cpp_Cbru,Cpp_Crub, Cbig_Cmol, Cbig_Cbru, Cbig_Crub, Cmol_Cbru, Cmol_Crub,Cbru_Crub)

library(ggplot2)
p_Dxy_nuc_all_1 <- ggplot(Dxy_nucleotide_all_1, aes(x=Dxy_nucleotide, fill=Name, colour=Name)) + geom_density(alpha=0.3)
p_Dxy_nuc_all_1

p_Dxy_nuc_all_1_colors <- p_Dxy_nuc_all_1 +scale_color_manual(values = c("Cpe_Cpp"="black", "Cbig_Cmol"="cyan", "Cbig_Cbru"="cyan1","Cbig_Crub"="cyan2", "Cmol_Cbru"="cyan3", "Cmol_Crub"="cyan4","Cbru_Crub"="deepskyblue","Cpp_Cbig"="chocolate", "Cpe_Cbig"="chocolate1", "Cpp_Cmol"="chocolate2","Cpe_Cmol"="chocolate3","Cpp_Cbru"="chocolate4","Cpe_Cbru"="coral","Cpp_Crub"="coral1","Cpe_Crub"="coral4"))+
  scale_fill_manual(values = c("Cpe_Cpp"="black", "Cbig_Cmol"="cyan", "Cbig_Cbru"= "cyan1" ,"Cbig_Crub"="cyan2", "Cmol_Cbru"="cyan3", "Cmol_Crub"="cyan4","Cbru_Crub"="deepskyblue","Cpp_Cbig"="chocolate", "Cpe_Cbig"="chocolate1", "Cpp_Cmol"="chocolate2","Cpe_Cmol"="chocolate3","Cpp_Cbru"="chocolate4","Cpe_Cbru"="coral","Cpp_Crub"="coral1","Cpe_Crub"="coral4"))
p_Dxy_nuc_all_1_colors

p_Dxy_nuc_all_1_zoomed<- p_Dxy_nuc_all_1_colors+coord_cartesian(xlim = c(0,0.01))
p_Dxy_nuc_all_1_zoomed

# Histogram
range<-max(Dxy_nucleotide_all_1$Dxy_nucleotide)-min(Dxy_nucleotide_all_1$Dxy_nucleotide)
p_hist_count_1<-ggplot(Dxy_nucleotide_all_1,aes(x=Dxy_nucleotide))+geom_histogram(fill="white",colour="black", binwidth = range/300)
p_hist_count_1
p_hist_count_1_zoomed<-p_hist_count_1+coord_cartesian(xlim = c(0,0.01))
library(ggpubr)
compare<- ggarrange(p_Dxy_nuc_all_1_colors, p_hist_count_1, ncol = 2, common.legend = FALSE,legend = FALSE, labels = "AUTO", font.label = list(size=18, face="plain"))
compare
compare_2<- ggarrange(p_Dxy_nuc_all_1_zoomed, p_hist_count_1_zoomed, ncol = 2, common.legend = FALSE,legend = FALSE, labels = "AUTO", font.label = list(size=18, face="plain"))
compare_2
p_hist_dens_1<-ggplot(Dxy_nucleotide_all_1,aes(x=Dxy_nucleotide))+geom_histogram(aes(y=..density..),
                                                                    fill="white",colour="black", binwidth = range/300)
p_hist_dens_1
p_hist_dens_1+coord_cartesian(xlim = c(0,0.01))


### Haplotype


# Cpe_Cpp_hap<- hapdiv_between_1[,1, drop=FALSE]
# colnames(Cpe_Cpp_hap) <- "Dxy_haplotype"
# Cpe_Cpp_hap$Name <- "Cpe_Cpp"
#
# Cpe_Cbig_hap<- hapdiv_between_1[,2, drop=FALSE]
# colnames(Cpe_Cbig_hap) <- "Dxy_haplotype"
# Cpe_Cbig_hap$Name <- "Cpe_Cbig"
#
# Cpe_Cmol_hap<- hapdiv_between_1[,3, drop=FALSE]
# colnames(Cpe_Cmol_hap) <- "Dxy_haplotype"
# Cpe_Cmol_hap$Name <- "Cpe_Cmol"
#
# Cpe_Cbru_hap<- hapdiv_between_1[,4, drop=FALSE]
# colnames(Cpe_Cbru_hap) <- "Dxy_haplotype"
# Cpe_Cbru_hap$Name <- "Cpe_Cbru"
#
# Cpe_Crub_hap<- hapdiv_between_1[,5, drop=FALSE]
# colnames(Cpe_Crub_hap) <- "Dxy_haplotype"
# Cpe_Crub_hap$Name <- "Cpe_Crub"
#
# Cpp_Cbig_hap<- hapdiv_between_1[,6, drop=FALSE]
# colnames(Cpp_Cbig_hap) <- "Dxy_haplotype"
# Cpp_Cbig_hap$Name <- "Cpp_Cbig"
#
# Cpp_Cmol_hap<- hapdiv_between_1[,7, drop=FALSE]
# colnames(Cpp_Cmol_hap) <- "Dxy_haplotype"
# Cpp_Cmol_hap$Name <- "Cpp_Cmol"
#
# Cpp_Cbru_hap<- hapdiv_between_1[,8, drop=FALSE]
# colnames(Cpp_Cbru_hap) <- "Dxy_haplotype"
# Cpp_Cbru_hap$Name <- "Cpp_Cbru"
#
# Cpp_Crub_hap<- hapdiv_between_1[,9, drop=FALSE]
# colnames(Cpp_Crub_hap) <- "Dxy_haplotype"
# Cpp_Crub_hap$Name <- "Cpp_Crub"
#
# Cbig_Cmol_hap<- hapdiv_between_1[,10, drop=FALSE]
# colnames(Cbig_Cmol_hap) <- "Dxy_haplotype"
# Cbig_Cmol_hap$Name <- "Cbig_Cmol"
#
# Cbig_Cbru_hap<- hapdiv_between_1[,11, drop=FALSE]
# colnames(Cbig_Cbru_hap) <- "Dxy_haplotype"
# Cbig_Cbru_hap$Name <- "Cbig_Cbru"
#
# Cbig_Crub_hap<- hapdiv_between_1[,12, drop=FALSE]
# colnames(Cbig_Crub_hap) <- "Dxy_haplotype"
# Cbig_Crub_hap$Name <- "Cbig_Crub"
#
# Cmol_Cbru_hap<- hapdiv_between_1[,13, drop=FALSE]
# colnames(Cmol_Cbru_hap) <- "Dxy_haplotype"
# Cmol_Cbru_hap$Name <- "Cmol_Cbru"
#
# Cmol_Crub_hap<- hapdiv_between_1[,14, drop=FALSE]
# colnames(Cmol_Crub_hap) <- "Dxy_haplotype"
# Cmol_Crub_hap$Name <- "Cmol_Crub"
#
# Cbru_Crub_hap<- hapdiv_between_1[,15, drop=FALSE]
# colnames(Cbru_Crub_hap) <- "Dxy_haplotype"
# Cbru_Crub_hap$Name <- "Cbru_Crub"
#
#
# Dxy_haplotype_all_1 <- rbind(Cpe_Cpp_hap, Cpe_Cbig_hap,Cpe_Cmol_hap, Cpe_Cbru_hap, Cpe_Crub_hap,Cpp_Cbig_hap,Cpp_Cmol_hap, Cpp_Cbru_hap,Cpp_Crub_hap, Cbig_Cmol_hap, Cbig_Cbru_hap, Cbig_Crub_hap, Cmol_Cbru_hap, Cmol_Crub_hap,Cbru_Crub_hap)
#
#
# p_Dxy_hap_all_1 <- ggplot(Dxy_haplotype_all_1, aes(x=Dxy_haplotype, fill=Name, colour=Name)) + geom_density(alpha=0.3)
# p_Dxy_hap_all_1<- p_Dxy_hap_all_1 +scale_color_manual(values = c("Cpe_Cpp"="black", "Cbig_Cmol"="cyan", "Cbig_Cbru"="cyan1","Cbig_Crub"="cyan2", "Cmol_Cbru"="cyan3", "Cmol_Crub"="cyan4","Cbru_Crub"="deepskyblue","Cpp_Cbig"="chocolate", "Cpe_Cbig"="chocolate1", "Cpp_Cmol"="chocolate2","Cpe_Cmol"="chocolate3","Cpp_Cbru"="chocolate4","Cpe_Cbru"="coral","Cpp_Crub"="coral1","Cpe_Crub"="coral4"))+
#   scale_fill_manual(values = c("Cpe_Cpp"="black", "Cbig_Cmol"="cyan", "Cbig_Cbru"= "cyan1" ,"Cbig_Crub"="cyan2", "Cmol_Cbru"="cyan3", "Cmol_Crub"="cyan4","Cbru_Crub"="deepskyblue","Cpp_Cbig"="chocolate", "Cpe_Cbig"="chocolate1", "Cpp_Cmol"="chocolate2","Cpe_Cmol"="chocolate3","Cpp_Cbru"="chocolate4","Cpe_Cbru"="coral","Cpp_Crub"="coral1","Cpe_Crub"="coral4"))
# p_Dxy_hap_all_1

### Second population set  pop1: Cpe-Cpp pop2: Cbig-Cmol-Cbru-Crub

Alignments_md10_popdefined_2 <- set.populations(Alignments_md10, list(c("Cpe_19esc","Cpe_21esc","Cpe_347bie","Cpe_348bie","Cpe_349bie","Cpp_25gab","Cpp_27gab", "Cpp_7aru", "Cpp_8aru", "Cpp_9aru"),
                                                                    c("Cbig_Bf1", "Cbig_Bf2", "Cbig_Bf3", "Cbig_Bm1", "Cbig_Bm2", "Cbig_Bm3", "Cbig_Bm4", "Cbig_Bm5","Cbig_Ef1","Cbig_Ef2","Cbig_Ef3","Cbig_Ef4","Cbig_Ef5", "Cbig_Em1","Cbig_Em2","Cbig_Em3","Cbig_Em4","Cbig_Em5",
                                                                    "Cmol_Bf1", "Cmol_Bf2", "Cmol_Bf3", "Cmol_Bf4", "Cmol_Bf5" ,"Cmol_Bm1", "Cmol_Bm2", "Cmol_Bm3", "Cmol_Bm4", "Cmol_Bm5","Cmol_Ef1","Cmol_Ef2","Cmol_Ef3","Cmol_Ef4","Cmol_Ef5", "Cmol_Em1","Cmol_Em2","Cmol_Em3","Cmol_Em4","Cmol_Em5",
                                                                    "Cbru_Bf1", "Cbru_Bf2", "Cbru_Bf3", "Cbru_Bf4", "Cbru_Bf5" ,"Cbru_Bm1", "Cbru_Bm2", "Cbru_Bm3", "Cbru_Bm4", "Cbru_Bm5","Cbru_Sf1", "Cbru_Sf2", "Cbru_Sf3", "Cbru_Sf4", "Cbru_Sf5" ,"Cbru_Sm1", "Cbru_Sm2", "Cbru_Sm3", "Cbru_Sm4", "Cbru_Sm5",
                                                                    "Crubratibialis_1", "Crubratibialis_2","Crubratibialis_3","Crubratibialis_4","Crubratibialis_5","Crubratibialis_6","Crubratibialis_7","Crubratibialis_8","Crubratibialis_9","Crubratibialis_10","Crubratibialis_12","Crubratibialis_14","Crubratibialis_15","Crubratibialis_16","Crubratibialis_19", "Crubratibialis_20")))
# Now calculate the Dxy again with these populations

Dxy_2<- diversity.stats.between(Alignments_md10_popdefined_2,subsites = FALSE, haplotype.mode = TRUE, nucleotide.mode = TRUE)

Dxy_2@populations # Defined populations beforehand: 1: Cpp+Cpe  2: Cbig+Cmol+Cbru+Crub

Dxy_2@nuc.diversity.between

nucdiv_between_2 <- as.data.frame(Dxy_2@nuc.diversity.between/Dxy_2@n.sites)
nucdiv_between_2
colMeans(nucdiv_between_2) ### 0.006258275
colnames(nucdiv_between_2) <- "Dxy_nuc"
nucdiv_between_2$Name <- "BIGvsPAR"

range<- max(nucdiv_between_2$Dxy_nuc)-min(nucdiv_between_2$Dxy_nuc)
p_hist_count<-ggplot(nucdiv_between_2,aes(x=Dxy_nuc))+geom_histogram(fill="white",colour="black", binwidth = range/300)
p_hist_count_zoomed<- p_hist_count+coord_cartesian(xlim = c(0,0.01))

p_hist_dens<-ggplot(nucdiv_between_2,aes(x=Dxy_nuc))+geom_histogram(aes(y=..density..),
                                     fill="white",colour="black", binwidth = range/300)

p_hist_dens+coord_cartesian(xlim = c(0,0.01))

BIGvsPAR <- ggplot(nucdiv_between_2, aes(x=Dxy_nuc, fill=Name, colour=Name)) + geom_density(alpha=0.3)+scale_fill_manual(values = c("BIGvsPAR"="black"))+scale_colour_manual(values = c("BIGvsPAR"="black"))
BIGvsPAR+geom_vline(xintercept = 0.006258275, lty=2)
BIGvsPAR_zoomed<- BIGvsPAR+coord_cartesian(xlim = c(0,0.01))

compare_3<- ggarrange(BIGvsPAR, p_hist_count, ncol = 2, common.legend = FALSE,legend = FALSE, labels = "AUTO", font.label = list(size=18, face="plain"))
compare_3
compare_4<- ggarrange(BIGvsPAR_zoomed, p_hist_count_zoomed, ncol = 2, common.legend = FALSE,legend = FALSE, labels = "AUTO", font.label = list(size=18, face="plain"))
compare_4
# write.table(nucdiv_between_2, file = "nucdiv_betweenBIGvsPAR", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
# library(readr)
# nucdiv_betweenBIGvsPAR <- read_delim("nucdiv_betweenBIGvsPAR",
#                               "\t", escape_double = FALSE, col_names = TRUE,
#                               trim_ws = TRUE)
