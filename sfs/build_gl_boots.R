#This script is used to combine gene SFS into a single 'genesum' SFS that contains all the data for a species comparison. It also builds 100 bootstrapped SFS using resampling with replacement of gene SFS.
#Run from the command line with Rscript command and the taxa pair abbreviation (cbig_b_cbru_b, etc.) as the first argument
#Alternatively, assign t manually in script.

options(scipen=999)

args = commandArgs(trailingOnly=TRUE)

t <- args[1]

  #create output directory
  dir.create(paste0('sfs/boot_sfs/',t), recursive = TRUE, showWarnings = FALSE)

  #Split taxa pair into taxa 1 and 2
  t1 <- substr(t,1,(nchar(t)/2)-0.5)
  t2 <- substr(t,(nchar(t)/2)+1.5,nchar(t))

  #read in bamlists to get number of individuals
  t1file <- read.csv(paste("lrz_",t1,".bamlist",sep = ""), header = FALSE)
  t2file <- read.csv(paste("lrz_",t2,".bamlist",sep = ""), header = FALSE)

  t1n <- length(t1file$V1)
  t2n <- length(t2file$V1)

  #read in filenames of gene SFS
  filelist <- list.files(paste0('sfs/per_gene_sfs/',t))

  #Sum gene SFS into total SFS and clean up.
  all_sfs <- t(sapply(filelist, function(x) scan(paste0('sfs/per_gene_sfs/',t,'/',x)), simplify = TRUE))

  all_sfs <- na.omit(all_sfs)

  summed_sfs <- colSums(all_sfs)

  if (t == 'cppar_cpery') {

    anc <- 'ancmol_'

  } else {

    anc <- ''

  }

  #Write combined gene SFS to single file with appendage 'genesum'
  sink(paste0('sfs/2dsfs_',t,'_p1_fold0_',anc,'genesum.sfs'), append = FALSE)
  cat(paste0((2*t1n)+1,' ',(2*t2n)+1,' unfolded "',t1,'" "',t2,'"'))
  cat('\n')
  sink()

  write(summed_sfs, file = paste0('sfs/2dsfs_',t,'_p1_fold0_',anc,'genesum.sfs'), ncolumns = length(summed_sfs), append = TRUE, sep = " ")

  #Build bootstrap SFS using sampling with replacement of gene SFS.
  for (b in 1:100) {

    sample_sfs <- all_sfs[sample(nrow(all_sfs),size=nrow(all_sfs),replace=TRUE),]
    boot_sfs <- colSums(sample_sfs)

    sink(paste0('sfs/boot_sfs/',t,'/boot',b,'.sfs'), append = FALSE)
    cat(paste0((2*t1n)+1,' ',(2*t2n)+1,' unfolded "',t1,'" "',t2,'"'))
    cat('\n')
    sink()

    write(boot_sfs, file = paste0('sfs/boot_sfs/',t,'/boot',b,'.sfs'), ncolumns = length(boot_sfs), append = TRUE, sep = " ")

  }
