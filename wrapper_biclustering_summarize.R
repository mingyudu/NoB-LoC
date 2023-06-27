
####### Create two summary file for protein and sample clusters, respectively

# Input:
### save.result: MCMC result
### burn.in: number of burn-in iterations
### protein_id: ID of proteins
### sample_id: ID of samples
### prot.clust.output: output file of protein cluster summary
### sample.clust.output: output file of sample cluster summary

# Output:
### infer.result: result of inference

bicluster_summarize <- function(save.result, burn.in, protein_id, sample_id, 
                                prot.clust.output, sample.clust.output){
  
  dyn.load("LC_pair.so")
  source('fn_LC_1_C_sig2_tau2.R')
  
  # burn.in <- 5000
  n.MCMC.iter <- dim(save.result$theta)[1]
  n.subject <- dim(save.result$theta)[2]
  n.protein <- dim(save.result$theta)[3]
  
  if (length(protein_id)!=n.protein){
    stop(cat('The length of protein_id should be equal to', n.protein))
  }
  if (length(sample_id)!=n.subject){
    stop(cat('The length of sample_id should be equal to', n.subject))
  }
  
  message('Calculate w.pair.prob ...')
  w.pair.prob <- fn.pair.prob.w(save.result$w[-(1:burn.in),], n.protein, (n.MCMC.iter - burn.in))
  
  best.w <- save.result$w[burn.in+w.pair.prob$sel,] # w_LS
  best.S <- max(best.w) # number of biclusters
  best.c <- array(NA, dim=c(best.S, n.subject)) # c_LS
  best.p <- best.K <- rep(NA, best.S)
  best.n <- array(0, dim=c(best.S, n.subject))
  
  message('Calculate c.pair.prob ...')
  for(i.s in 1:best.S)
  {
    gene.set <- (1:n.protein)[best.w==i.s]
    best.p[i.s] <- length(gene.set)
    c.pair.prob <- fn.pair.prob.c(gene.set, n.MCMC.iter, burn.in, n.subject, save.result$c, save.result$w)
    best.c[i.s, ] <- c.pair.prob$c
    best.K[i.s] <- max(c.pair.prob$c)
    
    for(i.k in 1:best.K[i.s])
    {
      best.n[i.s, i.k] <- sum(best.c[i.s, ]==i.k)
    }
    
  }
  best.m <- apply(best.n, 1, sum)
  
  ## summary of protein clusters
  
  message('Creating the summary for protein clusters...')
  num.clust <- max(best.w)
  sink(prot.clust.output)

  cat("Inactive protein clusters:\n")
  prots <- which(best.w==0)
  cat("Total number: ", length(prots), '\n')
  for (i in 1:length(prots)) {
    cat('\t', prots[i], '\t', protein_id[prots[i]], '\n')
  }

  for (i in 1:num.clust) {
    cat("Active protein bicluster ", i, ": \n")
    prots <- which(best.w==i)
    cat("Total number: ", length(prots), '\n')
    for (j in 1:length(prots)) {
      cat('\t', prots[j], '\t', protein_id[prots[j]], '\n')
    }
  }
  sink()

  ## summary of sample clusters

  message('Creating the summary for sample clusters...')
  sink(sample.clust.output)
  for (i in 1:num.clust) {
    cat("Bicluster ", i, ": \n")
    num.inact <- sum(best.c[i,]==0)
    if(num.inact==0){
      cat("\t", "**********There is no inactive sample.**********\n")
    }
    else{
      cat("\t", "# of inactive samples: ", num.inact, "\n")
      cat("\t", "Inactive samples: ", which(best.c[i,]==0), "\n")
    }
    num.sample.clust <- max(best.c[i,])
    cat("\t", "Total number of sample clusters: ", num.sample.clust, "\n")
    for (j in 1:num.sample.clust) {
      cat("\t", "Sample subcluster", j, ": \n")
      # cat("\t\t", which(best.c[i,]==j), "\n")
      num.sample <- which(best.c[i,]==j)
      for (k in num.sample) {
        cat("\t\t", k, '\t', sample_id[k], "\n")
      }
    }
  }
  sink()
  
  infer.result <- NULL
  infer.result$w.pair.prob <- w.pair.prob
  infer.result$best.w <- best.w
  infer.result$best.S <- best.S
  infer.result$best.c <- best.c
  infer.result$best.p <- best.p
  infer.result$best.K <- best.K
  infer.result$best.n <- best.n
  infer.result$best.m <- best.m
  return(infer.result)
}

save.result <- readRDS('./results/save.result.rds')
protein_id_list = readRDS('./data/protein_id_name.rds')
protein_id = protein_id_list$ID
sample_id = readRDS('./data/sample_id.rds')

infer.result = bicluster_summarize(save.result, burn.in = 5000, protein_id, sample_id, 
                    prot.clust.output = 'prot.txt', sample.clust.output = 'sample.txt')
saveRDS(infer.result, './results/infer.result.rds')
