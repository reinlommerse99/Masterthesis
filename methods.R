rm(list=ls())
library(PlackettLuce)
library(gosset)
library(climatrends)
library(chirps)
library(ggplot2)
library('plot.matrix')
library(bnlearn)
library(networktree)
library(partykit)
library(lubridate)
library(elevatr)


#Function to create data of rankings per farmer. Each row corresponds to a ranking of a particular farmer for a particular variety
#Also corresponds to the L in my thesis. See my thesis for a further explanation
L_maker <- function(covar, J = 3){
  L <- data.frame()
  
  #Extract all traits
  traits_worst <- colnames(covar)[grepl("_worst", colnames(covar))]
  traits_best <- colnames(covar)[grepl("_best", colnames(covar))]
  traits <- gsub("_worst", "", traits_worst)
  
  #Get the columns in which the varieties are mentioned
  varieties <- colnames(covar)[grepl("variety", colnames(covar))]
  
  tricot <- c("A", "B", "C")
  
  #Create data frame that contains numerical representation of rankings
  #A 1 corresponds to the most preferred, while a 3 corresponds to the inferior variety
  for (i in 1:nrow(covar)) {
    for (m in 1:length(traits)) {
      for (j in 1:length(varieties)) {
        
        #Give an id for the farmer. Best to do it yourself
        L[(J*(i-1) + j),1] <- i #covar$id[i]
        
        #State the variety
        L[(J*(i-1) + j),2] <- covar[,varieties[j]][i]
        
        #Collect the best and worst varieties
        best <- covar[,traits_best[m]][i]
        worst <- covar[,traits_worst[m]][i]
        
        #The scenario in which the farmer is indifferent
        if (is.na(best) & is.na(worst)) {
          L[(J*(i-1) + j),m+2] <- 0
          #The scenario in which only the worst variety is known
        } else if (is.na(best) & !is.na(worst)){
          if (worst == tricot[j]) {
            L[(J*(i-1) + j),m+2] <- 3
          } else{
            L[(J*(i-1) + j),m+2] <- 1.5
          }
          #The scenario in which only the best variety is known
        } else if (!is.na(best) & is.na(worst)){
          if (best == tricot[j]) {
            L[(J*(i-1) + j),m+2] <- 1
          } else{
            L[(J*(i-1) + j),m+2] <- 2.5
          }
          #Normal procedure otherwise. Give best rank 1, worst 3, middle 2
        } else if (best == tricot[j]) {
          L[(J*(i-1) + j),m+2] <- 1
        } else if (worst == tricot[j]) {
          L[(J*(i-1) + j),m+2] <- 3
        } else{
          L[(J*(i-1) + j),m+2] <- 2
        }
      }
    }
  }
  colnames(L) <- c("id", "Variety", traits)
  
  L <- L[rowSums(L[,-c(1,2)]) != 0,]
  
  L <- L[order(L$id),]
  return(L)
}


#Function that turns L into Y (thesis). Y contains of Plackett-Luce worth
Y_maker <- function(L, J=3){
  
  traits <- colnames(L[,-c(1,2)])
  L$id <- rep(1:(length(L$id)/J),each = J)
  
  #Create a list of rankings for Plackett-Luce (PL)
  R <- vector(mode = "list", length = length(traits))
  
  for (m in seq_along(traits)) {
    L_m <- L[c(1,2,(m+2))]
    
    R[[m]] <- rank_numeric(data = L_m,
                           items = "Variety",
                           input = traits[m],
                           id = "id",
                           ascending = TRUE)
  }
  
  #Apply PL to the rankings and create a data frame. Columns correspond with traits and rows with varieties
  mod <- list()
  
  #Get PL worth. Default value for method resulted in errors, for this reason we use BFGS.
  for (m in 1:length(R)) {
    mod[[m]] <- PlackettLuce(R[[m]], method = "BFGS")
  }
  
  Y <- data.frame()
  
  for (j in 1:length(mod)) {
    for (i in 1:length(mod[[j]]$coefficients)) {
      Y[i, j] <- mod[[j]]$coefficients[i]
    }
  }
  
  #As ties are prevalent in the data, we want to filter out the parameters corresponding to these ties
  variety_names <-c(row.names(as.data.frame(mod[[1]]$coefficients)),"tie2","tie3")
  
  Y <- Y[!(grepl("^tie*",variety_names)),]
  variety_names <- variety_names[!(grepl("^tie*",variety_names))]
  colnames(Y) <- traits
  
  Y$seeds <- variety_names
  
  return(Y)
}


#Create data set with all Kendall Tau coefficients
#This data will be used for the conditional inference trees
tau_maker <- function(L, J = 3){
  traits <- colnames(L[,-c(1,2)])
  taus <- data.frame()
  
  for (i in 1:(length(L$id)/J)) {
    n = 0
    for (j in 1:length(traits)) {
      for (m in 1:length(traits)) {
        if (m > j) {
          n = n + 1
          taus[i,n] <- cor(x = L[(J*(i-1) +1): (J*i), j + 2], y = L[(J*(i-1) +1): (J*i),m + 2], method = "kendall")
        }
      }
    }
  }
  taus[is.na(taus)] <- 0
  
  #Add names to dataset taus so we know which Kendall tau coefficient we see
  names <- c()
  for (j in 1:length(traits)) {
    for (m in 1:length(traits)) {
      if (m >j){
        names <- append(names, paste(traits[j], traits[m], sep = "_"))
      }
    }
  }
  
  colnames(taus) <- names
  taus$id <- L$id[seq(1,length(L$id),J)]
  return(taus)
}


X_maker <- function(L,J = 3){
  
  #Data frame with only traits data
  L_traits <- L[,-c(1,2)]
  traits <- colnames(L_traits)
  
  #Create data frame with +0.1 and -0.1 with the same dimensions as L_traits
  #This will be added to it set ties like the number 2.5 to either 2.4 or 2.6
  #This number will then be rounded to the nearest integer, making it either 2 or 3
  #This procedure mimics randomly assigning a variety that has either rank 2 or 3 to one of the options
  randommatrix <- matrix(sample(x = c(-0.1,0.1), size = nrow(L)*length(traits), replace = TRUE),
                         byrow = TRUE, ncol = length(traits))
  randommatrix <- data.frame(randommatrix)
  colnames(randommatrix) <- traits
  
  rounded_L_traits <- round(randommatrix + L_traits)
  
  #Index of the first variety per farmer
  f_idx <- seq(1,nrow(L),J)
  
  #This is added to randomly select a variety per farmer
  r_idx <- sample(x = seq(0, (J-1)), size = length(f_idx), replace = TRUE)
  
  #idxs are indices that select exactly one variety per farmer at random
  idxs <- f_idx + r_idx
  
  #The data in which every row contains the ranking of one variety per farmer  
  X <- rounded_L_traits[idxs,]
  
  #Assign a random ranking to the ones who are completely indifferent
  indifferences = X[X == 0]
  X[X == 0] <- sample(x = 1:J, size = length(indifferences), replace = T)
  
  #Make factors of data so it can be read as discrete data
  X[,] <- lapply(X, as.factor)
  X$id <- L$id[seq(1,length(L$id),J)]
  
  return(X)
}


#Functon for making the BN given input W
#Can change the score manually and give a blacklist
#This blacklist tells the BN which edges are deemed impossible
BN_maker <- function(W, score = "bic", bl = NULL){
  W <- W[,-ncol(W)]
  if (all(sapply(W, class) == "factor")) {
    dag <- tabu(W, score = score, blacklist = bl)
  } else if (all(sapply(W, class) == "numeric")) {
    dag <- tabu(W, score = "bge", blacklist = bl)
  } else{
    print("Not all vaiables are of the same class")
  }
  cp_dag <- cpdag(dag)
  return(dag)
}

#Both mle_CI and mle_ci_helper stem from a paper of Wiliams and Rast.
#Titled "Back to the basics: Rethinking partial correlation network methodology."
mle_CI <- function(X, alpha) {
  X <- as.matrix(X)
  if (!require("qgraph")){
    install.packages("qgraph")
  }
  if (!require("Matrix")) {
    install.packages("Matrix")
  }
  n <- nrow(X)
  p <- ncol(X)
  
  mle_cov <- crossprod(scale(X, scale = F)) / n
  
  mle_inv <- solve(mle_cov)
  
  par_cors <- as.matrix(qgraph::wi2net(mle_inv))
  mle_parcors <- mle_ci_helper(alpha = alpha, par_cors = par_cors, n = n, s = p-1)
  mle_inv <- mle_parcors$sig_mat * mle_inv
  list(mle_parcors = mle_parcors, mle_inv = mle_inv)
}

mle_ci_helper <- function(alpha, par_cors, s, n){
  mat <- matrix(0, nrow = s +1, ncol = s + 1)
  CI_ls <- list()
  par_cor <- par_cors[upper.tri(par_cors)]
  cov <- list()
  for (i in 1:length(par_cor)) {
    z_crit <- qnorm(1 - alpha/2)
    se <- sqrt(1/((n-2-3)))
    z <- log((1 + par_cor[i])/(1 - par_cor[i]))/2
    Z_L <- z - z_crit*se
    Z_U <- z + z_crit*se
    rho_L <- (exp(2*Z_L)-1)/(exp(2*Z_L)+1)
    rho_U <- (exp(2*Z_U)-1)/(exp(2*Z_U)+1)
    CI <- c(rho_L, rho_U)
    CI_ls[[i]] <- CI
    cov[[i]] <- ifelse(CI[1] < 0 & CI[2] >0, 0, 1)
  }
  
  ci_dat <- do.call(rbind.data.frame, CI_ls)
  colnames(ci_dat) <- c("low", "up")
  ci_dat$pcor <- unlist(par_cor)
  diag(mat) <- 1
  mat[upper.tri(mat)] <- unlist(cov)
  mat <- as.matrix(Matrix::forceSymmetric(mat))
  list(sig_mat = mat, par_cors = par_cors, par_sig = mat*par_cors, cis =ci_dat, cov_prob = unlist(cov))
}


#Get the partial correlation matrix if continuous data Y is the input
pcor_mat_makerY <- function(Y, alpha = 0.05){
  Y <- Y[,-ncol(Y)]
  traits <- colnames(Y)
  
  est_mle <- mle_CI(Y, alpha = alpha)
  
  #Obtain significant correlations
  mat <- est_mle$mle_parcors$par_sig
  return(mat)
}

#Get the partial correlation matrix if discrete data X is the input
#As X can be random, we bootstrap it to create a more robust result
pcor_mat_makerX <- function(L, alpha = 0.05, sim = 500, thresh_edge = 0.85){
  traits <- colnames(L)[-c(1,2)]
  
  mat_count <- matrix(rep(0,length(traits)*length(traits)),byrow = T,
                      ncol = length(traits))
  
  colnames(mat_count) <- rownames(mat_count) <- traits
  
  mat_par <- mat_count
  
  for (u in 1:sim) {
    X  <- X_maker(L)
    
    X <- X[,-ncol(X)]
    
    X[,] <- sapply(X, as.numeric)
    
    est_mle <- mle_CI(X, alpha = alpha)
    
    #Obtain significant correlations
    mat_count <- mat_count + est_mle$mle_parcors$sig_mat
    mat_par <- mat_par + est_mle$mle_parcors$par_cors
  }
  
  mat <- mat_par/sim
  mat[mat_count <= (thresh_edge*sim)] <- 0
  
  return(mat)
}

#Create a graph in which each edge weight is the Kendall Tau coefficient. 
kcor_mat_maker <- function(L, J=3, alpha = 0.05){
  traits <- colnames(L[,-c(1,2)])
  mat <- matrix(rep(1,(length(traits)*length(traits))), byrow = T, ncol = length(traits))
  colnames(mat) <- row.names(mat) <- traits
  
  for (m in 1:length(traits)) {
    for (n in 1:length(traits)) {
      if (n > m) {
        sum_taus = matrix(rep(0,(length(L$id)/J)))
        for (i in 1:(length(L$id)/J)) {
          add_tau <- cor(x = L[(J*(i-1) +1): (J*i), m + 2], y = L[(J*(i-1) +1): (J*i),n + 2], method = "kendall")
          if (!is.na(add_tau)) {
            sum_taus[i] <- add_tau
          }
        }
        
        avg_tau <- sum(sum_taus) / (length(L$id)/J)
        pval <- (sqrt((length(L$id)/J))*avg_tau)/sd(sum_taus)
        if (pval < qnorm((1 - alpha/2))) {
          mat[m,n] <- mat[n,m] <- 0
        } else{
          mat[m,n] <- mat[n,m] <- avg_tau
        }
      }
    }
  }
  return(mat)
}

#Create the CTRee-Rank. The input is the matrix of Kendall Tau coefficients for each farmer
ctree_maker <- function(taus, splitvars){
  N <- nrow(taus)
  
  taus <- taus[,-ncol(taus)]
  netdata <- taus; splitvars <- as.data.frame(splitvars)
  d <- cbind(netdata, splitvars)
  
  f1 <- as.formula(paste(c(paste(colnames(netdata),collapse=" + "), " ~ ", paste(colnames(splitvars), collapse=" + ")), collapse=""))
  
  ctree <- ctree(formula = f1, data = d, control = ctree_control(maxdepth = 2, minbucket = round(0.1*N)))
  return(ctree)
}

#Creates MOB trees. This is based on the paper of Jones
#Titled "Network trees: A method for recursively partitioning covariance structures."
mob_maker <- function(X, splitvars){
  X <- X[,-ncol(X)]
  X[,] <- sapply(X, as.numeric)
  mob <- networktree(nodevars=X,
                     splitvars=splitvars,
                     na.action=na.omit)
  return(mob)
}