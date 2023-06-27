########## MCMC implementation

####### Input:
### matrix: microbiome data: #sample by #protein
### pi_0: probability of being in any active protein set
### pi_1: probability of being an active sample given being in an active protein set
### beta: total mass parameter of protein clustering. The number of protein sets increases as beta increases
### M: total mass parameter of sample sub-clustering. The number of sample sets increases as M increases
### n.MCMC.iter: total number of MCMC iterations
### burn.in: number of burn-in iterations
### seed: random seed number

####### Output:
### save.result: MCMC result with n.MCMC.iter iterations

bicluster_mcmc <- function(matrix, pi_0, pi_1, beta, M, n.MCMC.iter, burn.in, seed){
  
  library(mnormt)
  source("fn_LC_1_C_sig2_tau2.R")
  dyn.load("LC_C_update_z_c_sig2.so")
  dyn.load("LC_C_update_w_sig2.so")
  dyn.load("LC_C_update_theta_sig2_tau2.so")
  
  if(!is.matrix(matrix)){
    matrix = as.matrix(matrix)
  }
  
  n.subject <- nrow(matrix); n.protein <- ncol(matrix)
  cat('Total number of samples is', n.subject, '.\n')
  cat('Total number of protein is', n.protein, '.\n')
  
  rownames(matrix)<-paste0('sample', 1:n.subject)
  colnames(matrix)<-paste0('protein', 1:n.protein)
  
  #Hierarchical clustering################
  distance_mat <- dist(t(matrix), method = 'euclidean')
  set.seed(240)
  Hierar_cl <- hclust(distance_mat, method = "complete")
  fit <- cutree(Hierar_cl, k = 30)
  tab <- table(fit)
  # table(fit)
  
  group.clust <- which(tab!=1)
  ini.w <- rep(0, n.protein)
  for (i in 1:length(group.clust)){
    ini.w[which(fit==group.clust[i])] <- i
  }
  ini.S <- length(group.clust)
  cat("Initial active protein clusters: ", ini.S, '\n')
  cat("Initial inactive proteins: ", sum(tab==1), '\n')
  
  #HYPER-PARAMETER VALUES SETTING######################################
  hyper <- NULL
  #GENE SETS
  ####PROBABILITY OF BEING IN ANY GENE SET.
  # hyper$pi_0 <- 0.01
  hyper$pi_0 <- pi_0
  # hyper$beta <- 1
  hyper$beta <- beta
  
  #DP FOR CLUSTERING GIVEN AN ACTIVE GENE SET
  # hyper$pi_1 <- 0.5
  hyper$pi_1 <- pi_1
  # hyper$M <- 1 
  hyper$M <- M
  
  #MEAN AND VAR FOR GENE IN ACTIVE GENE SETS AND SUBJECTS IN DP
  hyper$mu0 <- apply(matrix, 2, median)
  hyper$a0.g <- rep(NA, n.protein)
  hyper$b0.g <- rep(NA, n.protein)
  
  
  #MEAN AND VAR FOR GENE IN ACTIVE GENE SETS AND SUBJECTS NOT IN DP
  hyper$mu1 <- apply(matrix, 2, median)
  hyper$a1.g <- rep(NA, n.protein)
  hyper$b1.g <- rep(NA, n.protein)
  
  #MEAN AND VAR FOR GENE NOT IN ACTIVE GENE SETS AND ALL THE SUBJECTS
  hyper$mu2 <-  apply(matrix, 2, median)
  hyper$a2.g <- rep(NA, n.protein)
  hyper$b2.g <- rep(NA, n.protein)
  
  m <- apply(matrix, 2, var)
  v <- 1
  
  hyper$a0.g <- hyper$a1.g <- hyper$a2.g <- m^2/v + 2
  hyper$b0.g <- hyper$b1.g <- hyper$b2.g <- m*(m^2/v + 1)
  
  
  #HYPER-PARAMETERS FOR sig_g.
  hyper$a.g <- rep(NA, n.protein)
  hyper$b.g <- rep(NA, n.protein)
  
  # trim
  
  tmp <- matrix
  gupper.q <- apply(tmp, 2, quantile, probs=0.95)
  glower.q <-apply(tmp, 2, quantile, probs=0.05)
  
  v <- 1
  
  for(i.p in 1:n.protein)
  {
    tmp.exp <- matrix[,i.p]
    
    ind.1 <- (tmp.exp > gupper.q[i.p])
    ind.2 <- (tmp.exp < glower.q[i.p])
    
    m <- var(tmp.exp[!(ind.1)&!(ind.2)])
    
    hyper$a.g[i.p] <-  m^2/v + 2
    hyper$b.g[i.p] <- m*(m^2/v + 1)
    
  }
  
  ####ININITIALIZATION###########################
  set.seed(seed)
  cur.sample <- fn.initialization(hyper, matrix, ini.S, n.protein, n.subject, ini.w)
  
  #DO MCMC
  # n.MCMC.iter <- 35000
  # burn.in <- 5000
  n.MCMC.iter <- n.MCMC.iter
  burn.in <- burn.in
  
  save.result<- NULL
  save.result$theta <- array(NA, dim=c(n.MCMC.iter, n.subject, n.protein))
  save.result$w <- array(NA, dim=c(n.MCMC.iter, n.protein))
  save.result$S <- save.result$sig2 <- rep(NA, n.MCMC.iter)
  save.result$c <- save.result$n <- save.result$K <- save.result$m <- NULL
  save.result$sig2 <- save.result$tau20 <- save.result$tau21 <- save.result$tau22 <- array(NA, dim=c(n.MCMC.iter, n.protein))
  
  for(i.iter in 1:n.MCMC.iter)
  {
    #UPDATE Z AND C (SUBJECT CONFIGURATION FOR A GIVEN ACTIVE GENE SET)
    if(cur.sample$S > 0) {cur.sample <- fn.z.c.update.C(cur.sample, matrix, n.subject, hyper)}
    
    #UPDATE W (GENE CONFIGURATION)
    cur.sample <- fn.w.update.C(cur.sample, matrix, hyper, n.subject, n.protein)
    
    #UPDATE THETA AND SIG2
    tmp <- fn.theta.sig2.update.C(n.protein, n.subject, hyper, matrix, cur.sample)
    cur.sample$theta <- tmp$theta; cur.sample$sig2 <- tmp$sig2
    cur.sample$tau20 <- tmp$tau20; cur.sample$tau21 <- tmp$tau21; cur.sample$tau22 <- tmp$tau22
    
    save.result$theta[i.iter, , ] <- cur.sample$theta
    save.result$w[i.iter, ] <- cur.sample$w
    save.result$S[i.iter] <- cur.sample$S 
    save.result$sig2[i.iter,] <- cur.sample$sig2
    save.result$tau20[i.iter,] <- cur.sample$tau20
    save.result$tau21[i.iter,] <- cur.sample$tau21
    save.result$tau22[i.iter,] <- cur.sample$tau22
    
    if(cur.sample$S > 0)
    {
      save.result$c[[i.iter]] <- cur.sample$c
      save.result$n[[i.iter]] <- cur.sample$n
      save.result$K[[i.iter]] <- cur.sample$K
      save.result$m[[i.iter]] <- cur.sample$m
    }
    
    if((i.iter%%100)==0) print(paste("i.iter=", i.iter))
  }
  return(save.result)
}

library(readxl)
mat <- read_excel('./data/prot2peaks_data.xlsx')
mat <- subset(mat,select=-c(1:4))
mat <- t(as.matrix(mat))

n.subject <- nrow(mat); n.protein <- ncol(mat)
rownames(mat) <- paste0('subj', 1:n.subject)
colnames(mat) <- paste0('prot', 1:n.protein)

save.result = bicluster_mcmc(mat, pi_0 = 0.01, pi_1 = 0.5, beta = 1, M = 1, 
                     n.MCMC.iter = 10, burn.in = 5, seed = 345467)
saveRDS(save.result, './results/save.result.rds')
