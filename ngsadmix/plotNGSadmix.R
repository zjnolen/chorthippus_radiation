#This script was used to assemble plots for the structure analysis of the 84 individuals used in Zach's thesis.
#It is modified from the scripts built by Clara and Burcin.
#Redos of the analysis from Clara's thesis, which removed female Crub, can be found in the other file in this directory, plotNGSadmix_clara_redo.R


  colors = c("steelblue4","steelblue3","springgreen4","yellowgreen","slateblue4","thistle1","tomato3","orange1","orange4","yellow")
  mainlines = c(5,21,40,78)
  dotlines = c(1:88)
  spacevec=c(rep(0,10),1,rep(0,19),1,rep(0,17),1,rep(0,15),1,rep(0,19))


  #This chunk contains the color data for each analysis.
  #If the structure analysis is rerun, these numbers will change, as the cluster IDs are independent for each run.

    #for p1e-2, 50 reps
    colorscrambles[["50"]][["p1e-2"]][[2]] = c(3,1)
    colorscrambles[["50"]][["p1e-2"]][[3]] = c(1,5,3)
    colorscrambles[["50"]][["p1e-2"]][[4]] = c(5,1,7,3)
    colorscrambles[["50"]][["p1e-2"]][[5]] = c(5,1,3,8,7)
    colorscrambles[["50"]][["p1e-2"]][[6]] = c(2,5,8,1,7,3)
    colorscrambles[["50"]][["p1e-2"]][[7]] = c(1,2,9,5,7,3,8)
    colorscrambles[["50"]][["p1e-2"]][[8]] = c(4,3,8,1,2,5,9,7)
    colorscrambles[["50"]][["p1e-2"]][[9]] = c(6,3,1,4,7,5,9,2,8)
    colorscrambles[["50"]][["p1e-2"]][[10]] = c(2,6,3,7,1,4,10,5,9,8)

#The final chunk is the plotting chunk, the only thing that must be changed here is two variables: 'pval' and 'reps'.
#Other changes can be made as well to customize the plot output dimensions and name.

    pval = "p1e-2"
    reps = 50

    svg(file = paste("plots/NGSadmix",pval,reps,"reps.svg", sep = "_"), width = 9, height = 6)
    par(mar=c(0.75,0,0,0),oma=c(0,0,1,0))
    layout(c(1:10)) #10 is max clusters
    for (i in 2:10) { #10 is max clusters
      readin <- read.table(paste0("outputs_",reps,"reps/beagle_",pval,"_K",i,"_",reps,"reps.qopt"))
      reorder <- readin[c(80:84,75:79,39:58,1:18,59:74,29:38,19:28),]
      Q <- t(as.matrix(reorder))


      barplot(Q, col=colors[colorscrambles[[as.character(reps)]][[pval]][[i]]], border=colors[colorscrambles[[as.character(reps)]][[pval]][[i]]], axes=F, space=spacevec)
      sapply(mainlines, function(x) {lines(c(x,x),c(0,1),lty=1,lwd=1.5)})
      sapply(dotlines, function(x) {lines(c(x,x),c(0,1),lty=3,col="black")})
      mtext(paste0("K=",i), side=2, line=-2, font=2)
    }
    dev.off()
