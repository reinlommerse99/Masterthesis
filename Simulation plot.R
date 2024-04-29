## This function is useful for visualizing results
power_plot <- function(mob, ctree, ctree_rank, guide, ncheck=F, changeStat="tau"){
  Ns <- unique(guide$N)
  stepsizes <- unique(guide$stepsize)
  means_mob <- means_ctree <- means_ctree_rank <- n_sample_mob <- n_sample_ctree <- n_sample_ctree_rank <- vector()
  Ntemp <- rep(c(Ns[1], Ns[2], Ns[3]),3)
  stepTemp <- c(rep(stepsizes[1],3),rep(stepsizes[2],3),rep(stepsizes[3],3))
  sits <- c("sit1", "sit2", "sit3", "sit4", "sit5", "sit6",
            "sit7", "sit8", "sit9", "sit10")
  corChange <- c(0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  cols <- hcl.colors(3, "Geyser")   
  
  listdata <- list()
  for(i in 1:9){
    frametemp_mob <- mob[mob$N==Ntemp[i] & mob$stepsize==stepTemp[i],]
    frametemp_mob <- na.omit(frametemp_mob)
    frametemp_ctree <- ctree[ctree$N==Ntemp[i] & ctree$stepsize==stepTemp[i],]
    frametemp_ctree <- na.omit(frametemp_ctree)
    frametemp_ctree_rank <- ctree_rank[ctree_rank$N==Ntemp[i] & ctree_rank$stepsize==stepTemp[i],]
    frametemp_ctree_rank <- na.omit(frametemp_ctree_rank)
    if(dim(frametemp_mob)[1]!=0){
      for(j in 1:length(unique(sits))){
        means_mob[j] <- mean(frametemp_mob$split[frametemp_mob$sit == unique(sits)[j]])
        means_ctree[j] <- mean(frametemp_ctree$split[frametemp_ctree$sit == unique(sits)[j]])
        means_ctree_rank[j] <- mean(frametemp_ctree_rank$split[frametemp_ctree_rank$sit == unique(sits)[j]])
        n_sample_mob[j] <- sum(frametemp_mob$sit == unique(sits)[j])
        n_sample_ctree[j] <- sum(frametemp_ctree$sit == unique(sits)[j])
        n_sample_ctree_rank[j] <- sum(frametemp_ctree_rank$sit == unique(sits)[j])
      }
      listdata[[i]] <- list(means_mob,means_ctree, means_ctree_rank,
                            n_sample_mob,n_sample_ctree,n_sample_ctree_rank)  
    }
  }
  
  lineplot <- function(x1){
    lines(unique(corChange),listdata[[x1]][[1]],col="black")
    lines(unique(corChange),listdata[[x1]][[2]], lty=5,col="black")
    lines(unique(corChange),listdata[[x1]][[3]], lty=3,col="black")
  }
  
  uncertaintyplot <- function(x1){
    lower_mob <- qbinom(0.025, listdata[[x1]][[4]], listdata[[x1]][[1]]) / listdata[[x1]][[4]]
    upper_mob <- qbinom(0.975, listdata[[x1]][[4]], listdata[[x1]][[1]]) / listdata[[x1]][[4]]
    lower_ctree <- qbinom(0.025, listdata[[x1]][[5]], listdata[[x1]][[2]]) / listdata[[x1]][[5]]
    upper_ctree <- qbinom(0.975, listdata[[x1]][[5]], listdata[[x1]][[2]]) / listdata[[x1]][[5]]
    lower_ctree_rank <- qbinom(0.025, listdata[[x1]][[6]], listdata[[x1]][[3]]) / listdata[[x1]][[6]]
    upper_ctree_rank <- qbinom(0.975, listdata[[x1]][[6]], listdata[[x1]][[3]]) / listdata[[x1]][[6]]
    
    polygon(
      x = c(unique(corChange), rev(unique(corChange))),
      y = c(lower_mob, rev(upper_mob)),
      col = paste0(cols[1], "50"),
      border = NA
    )
    polygon(
      x = c(unique(corChange), rev(unique(corChange))),
      y = c(lower_ctree, rev(upper_ctree)),
      col = paste0(cols[2], "50"),
      border = NA
    )
    
    polygon(
      x = c(unique(corChange), rev(unique(corChange))),
      y = c(lower_ctree_rank, rev(upper_ctree_rank)),
      col = paste0(cols[3], "50"),
      border = NA
    )
  }
  
  dev.off()
  layout(matrix(c(10,10,10,1:9,11,11,11), byrow = TRUE, ncol = 3),
         heights = c(.2,.3,.3,.3,.1))
  par(oma = rep(4.5,4), mar = rep(.3,4))
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box(); axis(2, las = 1)
  uncertaintyplot(7); lineplot(7)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box()
  uncertaintyplot(8); lineplot(8)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box()
  uncertaintyplot(9); lineplot(9)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box(); axis(2, las = 1)
  uncertaintyplot(4); lineplot(4)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box()
  uncertaintyplot(5); lineplot(5)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box()
  uncertaintyplot(6); lineplot(6)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box(); axis(2, las = 1); axis(1,at=corChange, cex.axis = 0.85)
  uncertaintyplot(1); lineplot(1)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box(); axis(1,at=c(0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5), cex.axis = 0.85)
  uncertaintyplot(2); lineplot(2)
  
  plot(NULL, ylim=c(0,1), xlim=c(0.05,0.5), axes=FALSE); box(); axis(1,at=c(0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5), cex.axis = 0.85)
  uncertaintyplot(3); lineplot(3)

  mtext("Power",side=2, line=27,cex=0.7)
  mtext("Power",side=2, line=27, at=1.6,cex=0.7)
  mtext("Power",side=2, line=27, at=2.75,cex=0.7)
  
  
  plot(NULL,axes=F,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
  legend("center",
         legend = c("MOB", "CTree", "CTree-Rank"),
         lty=c("solid","longdash","dotted"),
         col=c(cols),
         pt.cex = 1,
         horiz = TRUE,
         bg = "grey"
  )
  
  
  if(changeStat == "tau"){
    mtext(expression(paste(Delta, tau)),las=1, side = 1, line = -1.25, outer = TRUE,cex=0.8)
    mtext(expression(paste(Delta, tau)), side = 1, line = -1.25, outer = TRUE,cex=0.8, at=0.17)
    mtext(expression(paste(Delta, tau)),las=1, side = 1, line = -1.25, outer = TRUE,cex=0.8, at = 0.83)
  }
  
  
  mtext(paste("N =", as.character(Ns[2])),  side=1, line = 0.5, outer = TRUE, cex=.85)
  mtext(paste("N =", as.character(Ns[1])),  side=1, line = 0.5, outer = TRUE, at = .17, cex=.85)
  mtext(paste("N =", as.character(Ns[3])), side=1, line = 0.5, outer = TRUE, at = .83, cex=.85)
  
  mtext(expression(paste(Delta, mu,"=0.5")), side = 4, las = 1, line = 0.25, outer = TRUE, at = .45, cex=.85)
  mtext(expression(paste(Delta, mu,"=0.2")),  side = 4, las = 1, line = 0.25, outer = TRUE, at = .18, cex=.85)
  mtext(expression(paste(Delta, mu,"=1")), side = 4, las = 1, line = 0.25, outer = TRUE, at = .72, cex=.85)
  
}

power_plot(mob_df, ctree_df, ctree_rank_df, guide)