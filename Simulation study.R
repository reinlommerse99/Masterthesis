rm(list=ls())
library(gosset)
library(PlackettLuce)
library(climatrends)
library(chirps)
library(ggplot2)
library('plot.matrix')
library(bnlearn)
library(networktree)
library(qgraph)
library(partykit)
library(mvtnorm)

#Function simulating the data process of tricot.
data_simulation <- function(perms, new_perms, N, ratio, M = 2, tricot = 3, J = 15){
  seeds <- seq(1,J)
  traits <- seq(1,M)
  dat <- data.frame()
  for (i in 1:N) {
    #The varieties farmer i receives
    varieties <- sample(seq(1,J), tricot)
    dat[(tricot*(i-1)+1):(tricot*(i-1)+tricot),1] <- i
    dat[(tricot*(i-1)+1):(tricot*(i-1)+tricot),2] <- varieties
    for (m in 1:M) {
      if (i < (ratio*N)) {
        dat[(tricot*(i-1)+1):(tricot*(i-1)+tricot),(2+m)] <- rmvnorm(1,perms[[m]][varieties])[,]
      }
      else{
        dat[(tricot*(i-1)+1):(tricot*(i-1)+tricot),(2+m)] <- rmvnorm(1,new_perms[[m]][varieties])[,]
      }
    }
  }
  data_df <- dat[,-c(1,2)]
  colnames(data_df) <- c("trait_1", "trait_2")
  return(data_df)
}

#Create data set with all Kendall Tau coefficients
#This data will be used for the conditional inference trees
ctree_rank_data <- function(data_df, N, M, tricot){
  taus <- data.frame()
  
  for (i in 1:N) {
    n = 0
    for (j in 1:M) {
      for (m in 1:M) {
        if (m > j) {
          n = n + 1
          taus[i,n] <- cor(x = data_df[(tricot*(i-1) +1): (tricot*i), j], y = data_df[(tricot*(i-1) +1): (tricot*i),m], method = "kendall")
        }
      }
    }
  }
  
  #Add names to dataset taus so we know which kendall tau coefficient we see
  names <- c()
  for (j in 1:(ncol(data_df))) {
    for (m in 1:(ncol(data_df))) {
      if (m >j){
        names <- append(names, paste(colnames(data_df)[j], colnames(data_df)[m], sep = "_"))
      }
    }
  }
  colnames(taus) <- names
  return(taus)
}

#This is in fact the same as making X. This is necessary for MOB and the normal ctree
mob_data <- function(data_df, N, M, tricot){
  dat_disc <- data_df
  for (i in 1:N) {
    for (m in 1:M) {
      dat_disc[(tricot*(i-1)+1):(tricot*(i-1)+tricot),m] <- rank(-data_df[(tricot*(i-1)+1):(tricot*(i-1)+tricot),m])
    }
  }
  
  #Index of the first variety per farmer
  f_idx <- seq(1,N*tricot,tricot)
  
  #This is added to randomly select a variety per farmer
  r_idx <- sample(x = seq(0, (tricot-1)), size = length(f_idx), replace = TRUE)
  
  #idxs are indices that select exactly one variety per farmer at random
  idxs <- f_idx + r_idx
  
  #The data in which every row contains the ranking of one variety per farmer  
  data_disc_idx <- dat_disc[idxs,]
  return(data_disc_idx)
}

#Make the tree
tree_maker <- function(N,tree_dat, sit, tau_change, Nnoise = 0, ratio = 0.5, method = "ctree_rank", fake = FALSE){
  
  num_zeros <- round(N * ratio)
  num_ones <- N - num_zeros
  splitvars <- as.data.frame(c(rep(0, num_zeros), rep(1, num_ones)))
  
  colnames(splitvars) <- "Splitter"
  splitvars$Splitter <- as.factor(splitvars$Splitter)
  
  netdata <- as.data.frame(tree_dat)
  
  if (Nnoise > 1){  
    noise_mat <- matrix(rnorm(Nnoise*N), N, Nnoise)
    colnames(noise_mat) <- noise <- paste("noise", 1:Nnoise, sep="")
    splitvars <- cbind(splitvars, noise_mat)
  }
  
  d <- cbind(netdata, splitvars)
  f1 <- as.formula(paste(c(paste(colnames(netdata),collapse=" + "), " ~ ", paste(colnames(splitvars), collapse=" + ")), collapse=""))
  
  if (method == "ctree_rank") {
    tr <- ctree(formula = f1, data = d)
    split <- if(is.null(tr)){NA}else if(length(tr)>1){TRUE}else{FALSE}
  } else if (method == "mob") {
    tr <- networktree(nodevars=netdata,
                      splitvars=splitvars,
                      method = "mob")
    split <- if(length(tr)>1){TRUE}else{FALSE}
  } else if (method == "ctree") {
    tr <- networktree(nodevars=netdata,
                      splitvars=splitvars,
                      method = "ctree")
    split <- if(length(tr)>1){TRUE}else{FALSE}
  }
  
  res <- list(d = d, tau_change = tau_change, sit = sit, tr = tr, split = split, method = method,
              ratio = ratio, Nnoise = Nnoise)
  return(res)
}

#Test whether the trees for a simulation detects a split for given parameters
#Very repetitive piece of code in which each block does the same for a different
#interval of Kendall tau
testsimulation <- function(N, ratio,Nnoise, sim, stepsize, tricot=3, M = 2, J = 15){
  
  repetition <- 0
  
  res_list <- list()
  
  mus <- cumsum(rep(stepsize,J))
  
  #Set levels of Kendall Tau bins and counters for each bin
  tau1_low <- 0.025
  tau1_up <- 0.075
  count1 <- 0
  
  tau2_low <- 0.075
  tau2_up <- 0.125
  count2 <- 0
  
  tau3_low <- 0.125
  tau3_up <- 0.175
  count3 <- 0
  
  tau4_low <- 0.175
  tau4_up <- 0.225
  count4 <- 0
  
  tau5_low <- 0.225
  tau5_up <- 0.275
  count5 <- 0
  
  tau6_low <- 0.275
  tau6_up <- 0.325
  count6 <- 0
  
  tau7_low <- 0.325
  tau7_up <- 0.375
  count7 <- 0
  
  tau8_low <- 0.375
  tau8_up <- 0.425
  count8 <- 0
  
  tau9_low <- 0.425
  tau9_up <- 0.475
  count9 <- 0
  
  tau10_low <- 0.475
  tau10_up <- 0.525
  count10 <- 0
  
  while ((count1 < sim) | (count2 < sim) | (count3 < sim) | (count4 < sim)
         | (count5 < sim) | (count6 < sim) | (count7 < sim) | (count8 < sim)
         | (count9 < sim)| (count10 < sim)) {
    
    perms <- list()
    perms[[1]] <- mus
    perms[[2]] <- sample(mus)
    
    new_perms <- list() 
    new_perms[[1]] <- mus
    new_perms[[2]] <- sample(mus)
    
    #Calculate the kendall tau before and after change in permutation
    tau_before <- kendallTau(perms[[2]], mus)$kendallTau
    tau_after <- kendallTau(new_perms[[2]], mus)$kendallTau
    tau_change <- abs(tau_after - tau_before)
    
    if ((tau_change > tau1_low) & (tau_change < tau1_up)) {
      if (count1 >= sim) {
        next
      } else{
        count1 <- count1 + 1
        sit <- "sit1"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N,tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N,tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob 
      }
    } else if ((tau_change > tau2_low) & (tau_change < tau2_up)) {
      if (count2 >= sim) {
        next
      } else{
        count2 <- count2 + 1
        sit <- "sit2"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau3_low) & (tau_change < tau3_up)) {
      if (count3 >= sim) {
        next
      } else{
        count3 <- count3 + 1
        sit <- "sit3"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau4_low) & (tau_change < tau4_up)) {
      if (count4 >= sim) {
        next
      } else{
        count4 <- count4 + 1
        sit <- "sit4"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    }
    else if ((tau_change > tau5_low) & (tau_change < tau5_up)) {
      if (count5 >= sim) {
        next
      } else{
        count5 <- count5 + 1
        sit <- "sit5"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau6_low) & (tau_change < tau6_up)) {
      if (count6 >= sim) {
        next
      } else{
        count6 <- count6 + 1
        sit <- "sit6"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    }else if ((tau_change > tau7_low) & (tau_change < tau7_up)) {
      if (count7 >= sim) {
        next
      } else{
        count7 <- count7 + 1
        sit <- "sit7"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau8_low) & (tau_change < tau8_up)) {
      if (count8 >= sim) {
        next
      } else{
        count8 <- count8 + 1
        sit <- "sit8"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau9_low) & (tau_change < tau9_up)) {
      if (count9 >= sim) {
        next
      } else{
        count9 <- count9 + 1
        sit <- "sit9"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    } else if ((tau_change > tau10_low) & (tau_change < tau10_up)) {
      if (count10 >= sim) {
        next
      } else{
        count10 <- count10 + 1
        sit <- "sit10"
        repetition <- repetition + 1
        
        data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
        taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
        data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
        res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
        res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "ctree")
        res_mob <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = tau_change, ratio = ratio, Nnoise = Nnoise, method = "mob")
        
        res_list[[((repetition - 1)*3 + 1)]] <- res_ctree_rank
        res_list[[((repetition - 1)*3 + 2)]] <- res_ctree
        res_list[[((repetition - 1)*3 + 3)]] <- res_mob
      }
    }
    
    
  }
  return(res_list)
}

#Also test when there should not be a split.
testnochange <- function(N,ratio,Nnoise, sim, stepsize, tricot=3, M = 2){
  
  res_list <- list()
  sit <- "nosit"
  mus <- cumsum(rep(stepsize,J))
  
  for (i in 1:sim) {
    
    perms <- list()
    perms[[1]] <- mus
    perms[[2]] <- sample(mus)
    
    new_perms <- perms
    
    data_df <- data_simulation(perms = perms, new_perms = new_perms, N = N, M=M, tricot = tricot, J=J, ratio = ratio)
    taus <- ctree_rank_data(data_df = data_df, N = N, M=M, tricot = tricot)
    data_disc_idx <- mob_data(data_df = data_df, N=N, M=M, tricot = tricot)
    res_ctree_rank <- tree_maker(N = N,tree_dat = taus, sit = sit, tau_change = 0, ratio = ratio, Nnoise = Nnoise, method = "ctree_rank")
    res_ctree <- tree_maker(N = N, tree_dat = data_disc_idx, sit = sit, tau_change = 0, ratio = ratio, Nnoise = Nnoise, method = "ctree")
    res_mob <- tree_maker(N =N, tree_dat = data_disc_idx, sit = sit, tau_change = 0, ratio = ratio, Nnoise = Nnoise, method = "mob")
    
    res_list[[((i - 1)*3 + 1)]] <- res_ctree_rank
    res_list[[((i - 1)*3 + 2)]] <- res_ctree
    res_list[[((i - 1)*3 + 3)]] <- res_mob
  }
  
  return(res_list)
}

#Number of seeds, varieties and farmers and amount of simulations
tricot <- 3
J <- 15
M <- 2
N <- 1500
sim <- 100

#Number of false splitting variables
Nnoise <- 0

#Ratio between the one group and the other
ratio <- 0.5
stepsize <- 1

#Create grid for all testing across different values of N, ratio and stepsize
guide<-expand.grid(N=c(160,600,1200),
                   Nnoise=c(0,5),
                   stepsize=c(0.2,0.5,1),
                   ratio = c(0.3,0.5))

simulation_df <- data.frame()

for (k in 1:nrow(guide)) {
  print(k)
  res_list <- testsimulation(N = guide$N[k], ratio = guide$ratio[k],
                             Nnoise = guide$Nnoise[k], sim = sim,
                             stepsize = guide$stepsize[k])
  for (i in 1:length(res_list)) {
    len <- length(res_list)
    simulation_df[((k-1)*len + i),1] <- res_list[[i]]$split
    simulation_df[((k-1)*len + i),2] <- res_list[[i]]$sit
    simulation_df[((k-1)*len + i),3] <- res_list[[i]]$ratio
    simulation_df[((k-1)*len + i),4] <- res_list[[i]]$Nnoise
    simulation_df[((k-1)*len + i),5] <- res_list[[i]]$method
    simulation_df[((k-1)*len + i),6] <- guide$N[k]
    simulation_df[((k-1)*len + i),7] <- guide$stepsize[k]
  }
}

colnames(simulation_df) <- c("split", "sit", "ratio", "Nnoise", "method", "N", "stepsize")

saveRDS(simulation_df, file="simulationdf.RData")

ctree_rank_df <- simulation_df[simulation_df$method == "ctree_rank" & simulation_df$ratio == 0.3 & simulation_df$Nnoise == 0,]
ctree_df <- simulation_df[simulation_df$method == "ctree" & simulation_df$ratio == 0.3 & simulation_df$Nnoise == 0,]
mob_df <- simulation_df[simulation_df$method == "mob" & simulation_df$ratio == 0.3 & simulation_df$Nnoise == 0,]


guide_fake<-expand.grid(N=c(160,600,1200),
                        stepsize=c(0.2,0.5,1),
                        ratio = c(0.3,0.5))

sim <- 300
simulation_fake_df <- data.frame()

for (k in 1:nrow(guide_fake)) {
  print(k)
  res_list_fake <- testnochange(N = guide_fake$N[k], ratio = guide_fake$ratio[k],
                                sim = sim, stepsize = guide_fake$stepsize[k], Nnoise = 0)
  for (i in 1:length(res_list_fake)) {
    len <- length(res_list_fake)
    simulation_fake_df[((k-1)*len + i),1] <- res_list_fake[[i]]$split
    simulation_fake_df[((k-1)*len + i),2] <- res_list_fake[[i]]$sit
    simulation_fake_df[((k-1)*len + i),3] <- res_list_fake[[i]]$ratio
    simulation_fake_df[((k-1)*len + i),4] <- res_list_fake[[i]]$Nnoise
    simulation_fake_df[((k-1)*len + i),5] <- res_list_fake[[i]]$method
    simulation_fake_df[((k-1)*len + i),6] <- guide_fake$N[k]
    simulation_fake_df[((k-1)*len + i),7] <- guide_fake$stepsize[k]
  }
}
colnames(simulation_fake_df) <- c("split", "sit", "ratio", "Nnoise", "method", "N", "stepsize")

saveRDS(simulation_fake_df, file="simulationfakedf.RData")


ctree_rank_fake_df <- simulation_fake_df[simulation_fake_df$method == "ctree_rank",]
ctree_fake_df <- simulation_fake_df[simulation_fake_df$method == "ctree",]
mob_fake_df <- simulation_fake_df[simulation_fake_df$method == "mob",]

error_mob_df <- data.frame()
error_ctree_df <- data.frame()
error_ctree_rank_df <- data.frame()


for (k in 1:nrow(guide_fake)) {
  error_mob_df[k,1:3] <- guide_fake[k,]
  error_mob_df[k,4] <- mean(mob_fake_df$split[mob_fake_df$ratio == guide_fake$ratio[k] & mob_fake_df$N == guide_fake$N[k] & mob_fake_df$stepsize == guide_fake$stepsize[k]])
  error_mob_df[k,5] <- "mob"
  
  error_ctree_df[k,1:3] <- guide_fake[k,]
  error_ctree_df[k,4] <- mean(ctree_fake_df$split[ctree_fake_df$ratio == guide_fake$ratio[k] & ctree_fake_df$N == guide_fake$N[k] & ctree_fake_df$stepsize == guide_fake$stepsize[k]])
  error_ctree_df[k,5] <- "ctree"
  
  error_ctree_rank_df[k,1:3] <- guide_fake[k,]
  error_ctree_rank_df[k,4] <- mean(ctree_rank_fake_df$split[ctree_rank_fake_df$ratio == guide_fake$ratio[k] & ctree_rank_fake_df$N == guide_fake$N[k] & ctree_rank_fake_df$stepsize == guide_fake$stepsize[k]])
  error_ctree_rank_df[k,5] <- "ctree_rank"
}


