#INTIALIZE THE SAMPLES EXCEPT GAMMA STAR
#para <- hyper; y <- RPPA; SS <- 5; N.gene <- n.protein; N.sub <- n.subject
#para <- hyper; y <- Y; SS <- ini.S; N.gene <- n.protein; N.sub <- n.subject
fn.initialization <- function(para, y, SS, N.gene, N.sub, ini.w)
{
	ini.sample <- NULL
	d.g <- (runif(N.gene) < para$pi_0)

	ini.sample$S <- SS  #number of gene sets

#	temp.w <- hclust(dist(apply(y[,d.g], 2, mean)), "ave")
#	ini.sample$w <- rep(NA, N.gene); ini.sample$w[d.g==F] <- 0  #configuration of genes to compose gene sets. (w==0)<=> (not in any gene gene set)
#	ini.sample$w[d.g==T] <- cutree(temp.w, k=ini.sample$S.p)  #w indicates which gene is in which gene set.

#	ini.sample$w <- rep(0, N.gene); ini.sample$w[1:5] <- 1; ini.sample$w[6:8] <- 2; ini.sample$w[9:12] <- 3

	ini.sample$w <- ini.w; 

#	g1 <- c(19, 28, 36, 9, 6, 59, 54, 3, 35); g2 <- c(32, 31, 30, 7, 8, 16, 18, 37, 40, 48, 63, 10, 24, 53)
#	g3 <- c(2, 15, 58, 51, 25, 62); g4 <- c(1, 50, 55, 45, 41, 47, 27, 43, 61); g5 <- c(46, 57, 42, 14, 39, 12, 11, 20, 26, 60, 33, 34, 38, 13, 17)
#	ini.sample$w <- rep(0, N.gene); ini.sample$w[g1] <- 1; ini.sample$w[g2] <- 2; ini.sample$w[g3] <- 3;
#	ini.sample$w[g4] <- 4; ini.sample$w[g5] <- 5;

	ini.sample$G.p <- sum(ini.sample$w !=0)  #number of genes in gene sets.

	ini.sample$p <- rep(NA, ini.sample$S)  #number of genes in each gene set.
	for(s in 1:ini.sample$S) ini.sample$p[s] <- sum(ini.sample$w==s)

	ini.sample$tau20 <- para$b0.g/(para$a0.g-1)
	ini.sample$tau21 <- para$b1.g/(para$a1.g-1)
	ini.sample$tau22 <- para$b2.g/(para$a2.g-1)
	ini.sample$sig2 <- para$b.g/(para$a.g-1) #0.05  #median(apply(y, 2, var))

	ini.sample$z <- ini.sample$c <- array(0, dim=c(N.sub, ini.sample$S))
	ini.sample$m <- rep(0, ini.sample$S)  #the number of subjects having z_is=1
	ini.sample$K <- rep(0, ini.sample$S)  #number of subject clusters for each active gene set.

	ini.sample$n <- array(0, dim=c(N.sub, ini.sample$S))  #locations of every subject cluster for evry gene get. #size of subject clusters for every gene set.  For any active gene set s, gam.star and n are a length-K_s vector but for each bookkeeping, I defined as a length-N.sub vector.

	for(s in (1:SS))
	{
		ind.g <- (ini.sample$w==s)

		#subset genes in gene set s
		tau20.s <- ini.sample$tau20[ind.g]; tau21.s <-  ini.sample$tau21[ind.g]
		mu0.s <-  para$mu0[ind.g]; mu1.s <-  para$mu1[ind.g]
		y.s <- as.matrix(y[ , ind.g]); sig2.s <- ini.sample$sig2[ind.g]

		tmp.in <- fn.impute.z.c.iter1(ini.sample$p[s], para$M, sig2.s, mu0.s, tau20.s, para$pi_1, N.sub, y.s)
		for(j in 1:10) tmp.in <- fn.z.c.one.gs.C(mu0.s, mu1.s, tau20.s, tau21.s, sig2.s, para$pi_1, para$M, tmp.in$z, tmp.in$m, tmp.in$c, tmp.in$K, tmp.in$n, N.sub, ini.sample$p[s], y.s)

		ini.sample$z[,s] <- tmp.in$z
		ini.sample$c[,s] <- tmp.in$c
		ini.sample$m[s] <- tmp.in$m
		ini.sample$K[s] <- tmp.in$K
		ini.sample$n[,s] <- tmp.in$n
	}

	return(ini.sample)
}




#p.s <- ini.sample$p[s]; M <- para$M; sig2 <- sig2.s; mu0 <- mu0.s; tau20 <- tau20.s;y <- y.s; pi_1 <- para$pi_1

#INITIALIZE Z.S, C.S, M.S, K.S, N.S, GAM.S FOR ONE TIME
#START WITH SUBJECT 1, AND THEN GO OVER SUBJECTS SEQUENTIALLY.
#p.s <- 1; mu0 <- mu0.s; tau20 <- tau20.s; y <- y.s
fn.impute.z.c.iter1 <- function(p.s, M, sig2, mu0, tau20, pi_1, N.sub, y)
{
	z.s <- c.s <- n.s <- rep(0, N.sub)
	gam.s <- array(NA, dim=c(N.sub, p.s))

	#the first subject takes z_is=1 and makes the first cluster by itself
	z.s[1] <- 1
	m.s <- 1
	K.s <- 1
	n.s[1] <- 1
	c.s[1] <- 1
	
	c.is <- c.s[1]
	var.tmp <- 1/(1/tau20 + n.s[c.is]/sig2)
	mean.tmp <- var.tmp*(mu0/tau20 + y[c.s==c.is,]/sig2) 
	gam.s[c.is, ] <- rnorm(p.s, mean.tmp, sqrt(var.tmp))

	for(i.sub in (2:N.sub))
	{
		z.s[i.sub] <- ifelse(runif(1) < pi_1, 1, 0)

		if(z.s[i.sub]==1)
		{
			m.s <- m.s + 1

			prob.i <- rep(NA, (K.s+1))
			for(k in (1:K.s))
			{
				prob.i[k] <- log(n.s[k]/(M+m.s-1)) + dmnorm(y[i.sub,], gam.s[k,], diag(sig2, p.s), log=TRUE)
			}

			prob.i[K.s+1] <- log(M/(M+m.s-1)) + dmnorm(y[i.sub,], mu0, diag((sig2 + tau20), p.s), log=TRUE)

			prob.i <- prob.i - max(prob.i)
			c.is <- sample((1:(K.s+1)), 1, FALSE, prob= exp(prob.i))
			c.s[i.sub] <- c.is

			if(c.is > K.s)
			{
				#START A NEW SUBJECT CLUSTER
				n.s[c.is] <- 1
				K.s <- K.s + 1
				var.tmp <- 1/(1/tau20 + n.s[c.is]/sig2)
				mean.tmp <- var.tmp*(mu0/tau20 + y[c.s==c.is,]/sig2) 
				gam.s[c.is, ] <- rnorm(p.s, mean.tmp, sqrt(var.tmp))  #location for a new cluster
			}else{
				#JOIN ONE OF EXISTING CLUSTERS
				n.s[c.is] <- n.s[c.is] + 1
				y.sum <- y[c.s==c.is,]
				if(n.s[c.is] > 1) y.sum <- apply(as.matrix(y.sum), 2, sum)
				
				var.tmp <- 1/(1/tau20 + n.s[c.is]/sig2)
				mean.tmp <- var.tmp*(mu0/tau20 + y.sum/sig2) 
				gam.s[c.is,] <- rnorm(p.s, mean.tmp, sqrt(var.tmp)) #update the location for the cluster that subject i.sub joins.
			}
		}
	}

	return(list(z=z.s, m=m.s, K=K.s, n=n.s, c=c.s))
}



#N.gene <- n.protein; N.sub <- n.subject; para <- hyper; y <- RPPA; this.sample <- cur.sample
fn.theta.sig2.update.C <- function(N.gene, N.sub, para, y, this.sample)
{
#	void fn_theta_sig2_update(int *NN_gene, int *NN_sub, 
#						 double *mu0, double  *tau20, double *mu1, double *tau21, double *mu2, double *tau22,
#						 double *a_g, double *b_g,
#						 double *y, double *theta, double *sig2,
#						 int *w, int *n, int *c, int *K)  

	y <- array(y, dim=c(1, N.gene*N.sub))[1,]
	theta <- rep(0, N.gene*N.sub)
		
	if(this.sample$S > 0)
	{
		n <- array(this.sample$n, dim=c(1, N.sub*this.sample$S))[1,]
		cc <- array(this.sample$c, dim=c(1, N.sub*this.sample$S))[1,]
		output <- .C("fn_theta_sig2_update", NN_gene=as.integer(N.gene), NN_sub=as.integer(N.sub), mu0=as.double(para$mu0), tau20=as.double(this.sample$tau20), mu1=as.double(para$mu1), tau21=as.double(this.sample$tau21), mu2=as.double(para$mu2), tau22=as.double(this.sample$tau22), a_g=as.double(para$a.g), b_g=as.double(para$b.g), a0_g=as.double(para$a0.g), b0_g=as.double(para$b0.g), a1_g=as.double(para$a1.g), b1_g=as.double(para$b1.g), a2_g=as.double(para$a2.g), b2_g=as.double(para$b2.g), y=as.double(y), theta=as.double(theta), sig2=as.double(this.sample$sig2), w=as.integer(this.sample$w), n=as.integer(n), c=as.integer(cc), K=as.integer(this.sample$K)) 
	}else{
		output <- .C("fn_theta_sig2_update", NN_gene=as.integer(N.gene), NN_sub=as.integer(N.sub), mu0=as.double(para$mu0), tau20=as.double(this.sample$tau20), mu1=as.double(para$mu1), tau21=as.double(this.sample$tau21), mu2=as.double(para$mu2), tau22=as.double(this.sample$tau22), a_g=as.double(para$a.g), b_g=as.double(para$b.g), a0_g=as.double(para$a0.g), b0_g=as.double(para$b0.g), a1_g=as.double(para$a1.g), b1_g=as.double(para$b1.g), a2_g=as.double(para$a2.g), b2_g=as.double(para$b2.g), y=as.double(y), theta=as.double(theta), sig2=as.double(this.sample$sig2), w=as.integer(this.sample$w), n=as.integer(0), c=as.integer(0), K=as.integer(0)) 
	}
	
	
	theta <- array(output$theta, dim=c(N.sub, N.gene))
	return(list(theta=theta, sig2=output$sig2, tau20=output$tau20, tau21=output$tau21, tau22=output$tau22))	
}




#y.s <- as.matrix(y[, ind.g]); z.s <- this.sample$z[,i.gs]; m.s <- this.sample$m[i.gs]; c.s <- this.sample$c[,i.gs]; K.s <- this.sample$K[i.gs]; n.s <- this.sample$n[,i.gs];
#mu0.s <- hyper$mu0[ind.g]; tau20.s <- hyper$tau20[ind.g]; mu1.s <- hyper$mu1[ind.g]; tau21.s <- hyper$tau21[ind.g]; pi_1 <- hyper$pi_1
#M <- hyper$M; p.s <- this.sample$p[i.gs]; sig2.s <- this.sample$sig2[ind.g]

#sig2 <- ini.sample$sig2[ind.g]; M <- para$M; p.s <- ini.sample$p[s]
#z.s <- tmp.in$z.s; m.s <- tmp.in$m.s; c.s <- tmp.in$c.s; K.s <- tmp.in$K.s; n.s <- tmp.in$n.s
fn.z.c.one.gs.C <- function(mu0.s, mu1.s, tau20.s, tau21.s, sig2.s, pi_1, M, z.s, m.s, c.s, K.s, n.s, N.sub, p.s, y.s)
{

#	void fn_z_c_one_gs(double *mu_0, double *mu_1, double *tau2_0, double *tau2_1, double *sig2, double *ppi_1, double *MM, int *z, int *m, int *c, int *K, int *n, int *nn_sub, int*pp, double *y)
	y.s <- array(y.s, dim=c(1, N.sub*p.s))[1,]
	
	output <- .C("fn_z_c_one_gs", mu_0=as.double(mu0.s), mu_1=as.double(mu1.s), tau2_0=as.double(tau20.s), tau2_1=as.double(tau21.s), sig2=as.double(sig2.s), ppi_1=as.double(pi_1), MM=as.double(M), z=as.integer(z.s), m=as.integer(m.s), c=as.integer(c.s), K=as.integer(K.s), n=as.integer(n.s), nn_sub=as.integer(N.sub), pp=as.integer(p.s), y=as.double(y.s))
	
	return(list(z=output$z, m=output$m, c=output$c, K=output$K, n=output$n))
}



#this.sample <- cur.sample; y <- RPPA; N.sub <- n.subject
#this.sample <- cur.sample; y <- Y; N.sub <- n.subject
fn.z.c.update.C <- function(this.sample, y, N.sub, hyper)
{
	new.sample <- this.sample
	ind.gs <- 1
	
	for(i.gs in (1:this.sample$S))
#	for(i.gs in (1:3))
	{
		ind.g <- (this.sample$w==i.gs)  #subset genes in a active gene set s

		tmp <- fn.z.c.one.gs.C(hyper$mu0[ind.g], hyper$mu1[ind.g], this.sample$tau20[ind.g], this.sample$tau21[ind.g], this.sample$sig2[ind.g], hyper$pi_1, hyper$M, this.sample$z[,i.gs], this.sample$m[i.gs], this.sample$c[,i.gs], this.sample$K[i.gs], this.sample$n[,i.gs], N.sub, this.sample$p[i.gs], as.matrix(y[, ind.g]))

		new.sample$z[,i.gs] <- tmp$z  #size N.subject vector indicating which subject choses a gene set s
		new.sample$m[i.gs] <- tmp$m   #number of subjects who choose a gene set s
		new.sample$c[,i.gs] <- tmp$c  #size N.subject vetor indicating which subject is in which cluter.
		new.sample$K[i.gs] <- tmp$K   #number of cluters.
		new.sample$n[,i.gs] <- tmp$n   #size K vector representing cluster sizes.
	}
	
	return(new.sample)
	
}

#set.seed(43432)
#this.sample <- cur.sample; para <- hyper; y <- RPPA; N.sub <- n.subject; N.gene <- n.protein
#this.sample <- cur.sample; para <- hyper; y <- Y; N.sub <- n.subject; N.gene <- n.protein
fn.w.update.C <- function(this.sample, y, para, N.sub, N.gene)
{

#	void fn_w_update(double *pi_0, double *pi_1, double *Be, double *M,
#					 double *mu_0, double *tau2_0, double *mu_1, double *tau2_1, double *mu_2, double *tau2_2,
#					 double *sig2,
#					 int *w, int *p, int *G_p, int *S_p, int *S,
#					 int *z, int *n, int *c, int *K, int *m,
#					 double *y, int *nn_sub, int *nn_gene)

	y <- array(y, dim=c(1, N.sub*N.gene))[1,]
	
	#In R, I define z as a N.sub by S matrix
	#But in C, I define z as a N.sub by N.gene matrix which is its max dim.  
	#So I fill the rest of it as 0. the same for n and c
	if(this.sample$S > 0)
	{
		z <- array(this.sample$z, dim=c(1, N.sub*this.sample$S))[1,]; z <- c(z, rep(0, N.sub*(N.gene-this.sample$S)))
		n <- array(this.sample$n, dim=c(1, N.sub*this.sample$S))[1,]; n <- c(n, rep(0, N.sub*(N.gene-this.sample$S)))
		c <- array(this.sample$c, dim=c(1, N.sub*this.sample$S))[1,]; c <- c(c, rep(0, N.sub*(N.gene-this.sample$S)))

		K <- c(this.sample$K, rep(0, (N.gene - this.sample$S))); 
		m <- c(this.sample$m, rep(0, (N.gene - this.sample$S))) 
		p <- c(this.sample$p, rep(0, (N.gene - this.sample$S))) 

	}else{
		z <- rep(0, N.sub*N.gene)
		n <- rep(0, N.sub*N.gene)
		c <- rep(0, N.sub*N.gene)

		K <- rep(0, N.gene); 
		m <- rep(0, N.gene) 
		p <- rep(0, N.gene) 
	}


	output <- .C("fn_w_update", pi_0=as.double(para$pi_0), pi_1=as.double(para$pi_1), Be=as.double(para$beta), M=as.double(para$M), mu_0=as.double(para$mu0), tau2_0=as.double(this.sample$tau20), mu_1=as.double(para$mu1), tau2_1=as.double(this.sample$tau21), mu_2=as.double(para$mu2), tau2_2=as.double(this.sample$tau22), sig2=as.double(this.sample$sig2), w=as.integer(this.sample$w), p=as.integer(p), G_p=as.integer(this.sample$G.p), S=as.integer(this.sample$S), z=as.integer(z), n=as.integer(n), c=as.integer(c), K=as.integer(K), m=as.integer(m), y=as.double(y), nn_sub=as.integer(N.sub), nn_gene=as.integer(N.gene))

	this.sample$w <- output$w; this.sample$G.p <- output$G_p; this.sample$S <- output$S

	if(this.sample$S > 0)
	{
		this.sample$m <- output$m[1:output$S]; this.sample$K <- output$K[1:output$S];		this.sample$p <- output$p[1:output$S]; 

		z <- array(output$z, dim=c(N.sub, N.gene))[,1:this.sample$S]
		n <- array(output$n, dim=c(N.sub, N.gene))[,1:this.sample$S]
		c <- array(output$c, dim=c(N.sub, N.gene))[,1:this.sample$S]
		this.sample$z <- as.matrix(z); this.sample$n <- as.matrix(n); this.sample$c <- as.matrix(c)
	}else{
		this.sample$m <- NA; this.sample$K <- NA; this.sample$p <- NA
		this.sample$z <- NA; this.sample$n <- NA; this.sample$c <- NA
	}

	return(this.sample)
}


#SS <- save.result$w[2001:2003,]; N.trt <- n.protein; N.sample <- 3
fn.pair.prob.w <- function(SS, N.trt, N.sample)
{

	SS <- (array(t(SS), dim=c(1, N.trt*N.sample)))[1,]
	output1 <- .C("Calc_Pair_Prob", nn=as.integer(N.trt), TT=as.integer(N.sample), Ss_sample=as.integer(SS), Prob=as.integer(rep(0, N.trt*N.trt)))

	Prob <- output1$Prob/N.sample

	output2 <- .C("Calc_Squared_Distance", nn=as.integer(N.trt), TT=as.integer(N.sample), Ss_sample=as.integer(SS), Prob=as.double(Prob), sq_dis=as.double(rep(0, N.sample)))

	Prob <- array(Prob, dim=c(N.trt, N.trt))	
	sel <- which.min(output2$sq_dis)

	return(list(Prob=Prob, sel=sel, sq_dis=output2$sq_dis))
}


#gene.set <- gene.set; n.MCMC <- n.MCMC.iter; n.burn.in <- burn.in; N.sub <- n.subject; cc.all <- save.result$c; w <- save.result$w


fn.pair.prob.c <- function(gene.set, n.MCMC, n.burn.in, N.sub, cc.all, w)
{ 
	cc <- array(NA, dim=c((n.MCMC - n.burn.in), N.sub))
	t <- 1

	for(i.iter in (n.burn.in+1):n.MCMC)
	{
		w.tmp <- w[i.iter,]	
		
		s <- w.tmp[gene.set[1]]
		p.s <- length(gene.set)
	
		ind.1 <- (sum(w.tmp[gene.set]==s) == p.s) #all the genes in the same gene set 
		ind.2 <- (sum(w.tmp==s) == p.s)  #only the genes in the gene set
	
		if(ind.1 & ind.2)
		{
		   cc[t,] <- (cc.all[[i.iter]])[,s]
		   t <- t + 1
#		   print(i.iter)
		}
	}

	N.sample <- t-1; print(N.sample)
	cc <- cc[(1:N.sample),]
	cc <- (array(t(cc), dim=c(1, N.sub*N.sample)))[1,]

	output1 <- .C("Calc_Pair_Prob", nn=as.integer(N.sub), TT=as.integer(N.sample), Ss_sample=as.integer(cc), Prob=as.integer(rep(0, N.sub*N.sub)))

	Prob <- output1$Prob/N.sample

	output2 <- .C("Calc_Squared_Distance", nn=as.integer(N.sub), TT=as.integer(N.sample), Ss_sample=as.integer(cc), Prob=as.double(Prob), sq_dis=as.double(rep(0, N.sample)))

	Prob <- array(Prob, dim=c(N.sub, N.sub))	
	sel <- which.min(output2$sq_dis)
	sel <- N.sub*(sel - 1) + 1

	return(list(Prob=Prob, c=cc[(sel:(sel+N.sub-1))], sq_dis=output2$sq_dis))
	return(Prob)
}



