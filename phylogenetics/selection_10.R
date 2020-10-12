load("data-md10")
percentNs_md10$Cbig_Gf1 <- NULL
percentNs_md10$Cbig_Gf2 <- NULL
percentNs_md10<- percentNs_md10[,c(59,60,61,62,63,64,65,66,67,68,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
                                   22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84)]
## Extract the populations as separate data frames
Cpe_percentNs <- percentNs_md10[,c(1:5)]
Cpp_percentNs <- percentNs_md10[,c(6:10)]


Cbig_B_percentNs <- percentNs_md10[,c(11:18)]
Cbig_E_percentNs <- percentNs_md10[,c(19:28)]

Cbru_B_percentNs <- percentNs_md10[,c(29:38)]
Cbru_S_percentNs <- percentNs_md10[,c(39:48)]

Cmol_B_percentNs <- percentNs_md10[,c(49:58)]
Cmol_E_percentNs <- percentNs_md10[,c(59:68)]


Crub_percentNs <- percentNs_md10[,c(69:84)]

library(readr)
GeneIDs_ORF <- read_delim("GeneIDs_ORF",
                                "\t", escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE)
## Look at the percentages of the genes which are 100% and more than 50% missing
## in the individuals
fullmissingpercent_CbigB<- ((rowSums(Cbig_B_percentNs==100.00000000))/8)*100
halfmissingpercent_CbigB<- ((rowSums(Cbig_B_percentNs>50.0000000))/8)*100
missing_genes_Cbig_B<- cbind(GeneIDs_ORF, halfmissingpercent_CbigB)

fullmissingpercent_CbigE<- ((rowSums(Cbig_E_percentNs==100.00000000))/10)*100
halfmissingpercent_CbigE<- ((rowSums(Cbig_E_percentNs>50.0000000))/10)*100
missing_genes_Cbig_E<- cbind(GeneIDs_ORF, halfmissingpercent_CbigE)

fullmissingpercent_CbruB<- ((rowSums(Cbru_B_percentNs==100.00000000))/10)*100
halfmissingpercent_CbruB<- ((rowSums(Cbru_B_percentNs>50.0000000))/10)*100
missing_genes_Cbru_B<- cbind(GeneIDs_ORF, halfmissingpercent_CbruB)

fullmissingpercent_CbruS<- ((rowSums(Cbru_S_percentNs==100.00000000))/10)*100
halfmissingpercent_CbruS<- ((rowSums(Cbru_S_percentNs>50.0000000))/10)*100
missing_genes_Cbru_S<- cbind(GeneIDs_ORF, halfmissingpercent_CbruS)

fullmissingpercent_CmolB<- ((rowSums(Cmol_B_percentNs==100.00000000))/10)*100
halfmissingpercent_CmolB<- ((rowSums(Cmol_B_percentNs>50.0000000))/10)*100
missing_genes_Cmol_B<- cbind(GeneIDs_ORF, halfmissingpercent_CmolB)

fullmissingpercent_CmolE<- ((rowSums(Cmol_E_percentNs==100.00000000))/10)*100
halfmissingpercent_CmolE<- ((rowSums(Cmol_E_percentNs>50.0000000))/10)*100
missing_genes_Cmol_E<- cbind(GeneIDs_ORF, halfmissingpercent_CmolE)

fullmissingpercent_Cpp<- ((rowSums(Cpp_percentNs==100.00000000))/5)*100
halfmissingpercent_Cpp<- ((rowSums(Cpp_percentNs>50.0000000))/5)*100
missing_genes_Cpp<- cbind(GeneIDs_ORF, halfmissingpercent_Cpp)

fullmissingpercent_Cpe<- ((rowSums(Cpe_percentNs==100.00000000))/5)*100
halfmissingpercent_Cpe<- ((rowSums(Cpe_percentNs>50.0000000))/5)*100
missing_genes_Cpe<- cbind(GeneIDs_ORF, halfmissingpercent_Cpe)

fullmissingpercent_Crub<- ((rowSums(Crub_percentNs==100.00000000))/16)*100
halfmissingpercent_Crub<- ((rowSums(Crub_percentNs>50.0000000))/16)*100
missing_genes_Crub<- cbind(GeneIDs_ORF, halfmissingpercent_Crub)

## subset the data by filtering out the the genes which are missing in more than
## 70% percent of the individuals

ss_missing_genes_CbigB_30 <- subset(missing_genes_Cbig_B, halfmissingpercent_CbigB<=70)
ss_missing_genes_CbigE_30 <- subset(missing_genes_Cbig_E, halfmissingpercent_CbigE<=70)

ss_missing_genes_CbruB_30 <- subset(missing_genes_Cbru_B, halfmissingpercent_CbruB<=70)
ss_missing_genes_CbruS_30 <- subset(missing_genes_Cbru_S, halfmissingpercent_CbruS<=70)
ss_missing_genes_CmolB_30 <- subset(missing_genes_Cmol_B, halfmissingpercent_CmolB<=70)
ss_missing_genes_CmolE_30 <- subset(missing_genes_Cmol_E, halfmissingpercent_CmolE<=70)
ss_missing_genes_Cpp_30 <- subset(missing_genes_Cpp, halfmissingpercent_Cpp<=70)
ss_missing_genes_Cpe_30 <- subset(missing_genes_Cpe, halfmissingpercent_Cpe<=70)
ss_missing_genes_Crub_30 <- subset(missing_genes_Crub, halfmissingpercent_Crub<=70)

library(dplyr)
common_genes_1_30 <- inner_join(ss_missing_genes_CbigB_30, ss_missing_genes_CbigE_30, by="X1")
common_genes_2_30 <- inner_join(common_genes_1_30, ss_missing_genes_CbruB_30, by="X1")
common_genes_3_30 <- inner_join(common_genes_2_30, ss_missing_genes_CbruS_30, by="X1")
common_genes_4_30 <- inner_join(common_genes_3_30, ss_missing_genes_CmolB_30, by="X1")
common_genes_5_30 <- inner_join(common_genes_4_30, ss_missing_genes_CmolE_30, by="X1")
common_genes_6_30 <- inner_join(common_genes_5_30, ss_missing_genes_Cpp_30, by="X1")
common_genes_7_30 <- inner_join(common_genes_6_30, ss_missing_genes_Cpe_30, by="X1")
all_common_30 <- inner_join(common_genes_7_30, ss_missing_genes_Crub_30, by="X1")

geneIDs_filt_md10_30p <- all_common_30[,1]
#write.table(geneIDs_filt_md10_30p, file = "geneIDs_filt_md10_30p",quote=FALSE, sep="\t", row.names=FALSE, col.names = FALSE)


percentNs_md10_wnames <- cbind(GeneIDs_ORF, percentNs_md10)
common_with_stat_30p<- percentNs_md10_wnames[percentNs_md10_wnames$X1 %in% geneIDs_filt_md10_30p, ]

common_with_stat_30p<-common_with_stat_30p %>% dplyr::rename(Crubratibialis_1 = Crub_1, Crubratibialis_2 = Crub_2,
                                                             Crubratibialis_3 = Crub_3,Crubratibialis_4 = Crub_4,
                                                             Crubratibialis_5 = Crub_5,Crubratibialis_6 = Crub_6,
                                                             Crubratibialis_7 = Crub_7,Crubratibialis_8 = Crub_8,
                                                             Crubratibialis_9 = Crub_9,Crubratibialis_10 = Crub_10,
                                                             Crubratibialis_12 = Crub_12,Crubratibialis_14 = Crub_14,
                                                             Crubratibialis_15 = Crub_15,Crubratibialis_16 = Crub_16,
                                                             Crubratibialis_19 = Crub_19,Crubratibialis_20 = Crub_20)

## For each of these genes, how many individuals I am gonna filter out?
## Count the NA values in each row

common_with_stat_30p$ex_count <- rowSums(common_with_stat_30p[-1] > 50.00)
common_with_stat_30p[common_with_stat_30p>50] <- "NA"
common_with_stat_30p$X1 <- geneIDs_filt_md10_30p

none<- common_with_stat_30p[which(common_with_stat_30p$ex_count==0, arr.ind = TRUE),]
none_geneIDs<- none[,1, drop=FALSE]

#write.table(none_geneIDs, file = "none_geneIDs_md10", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

not_none<- common_with_stat_30p[which(common_with_stat_30p$ex_count!=0, arr.ind = TRUE),]
rest_geneIDs <- not_none[,1,drop=FALSE]
View(rest_geneIDs)

#write.table(rest_geneIDs, file = "rest_geneIDs_md10", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

## Take the individuals names which are not NA for each gene.

not_none_2<- not_none[,1:85]
names<- apply(not_none_2[-1], 1, function(i) paste(names(i[i !="NA" ]), collapse = "\t"))
names<- data.frame(names)


setwd("genes_individuals_10")

for(i in 1:nrow(names)){
  myfile<-paste0(rest_geneIDs[i,], "_", "IDs")
  write.table(names[i,1],myfile,sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
}
