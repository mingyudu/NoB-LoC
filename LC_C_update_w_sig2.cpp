#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>

extern "C"
{
//	#include "LC_C_update_z_c.h"

	//WHEN UPDATING W, IF A GENE IS ACTIVE, THEN ADJUST OHTERS ACCORDING TO REMOVAL OF THAT GENE. 
	void fn_remove_gene(int *G_p, int *p, int *w, int w_g, int *S, int *z, int *n, int *c, int *K, int *m, int n_gene, int n_sub)
	{
		int s, g, i;
		
		(*G_p) = (*G_p) - 1;    //decrease the number of genes in any gene set by 1
		p[w_g-1] = p[w_g-1] - 1;    //decrease the number of genes in the gene set to which gene g belongs.

//		printf("\n REMOVE!!!! w_g=%d, p_g=%d", w_g, p[w_g-1]);
		
		
		if(p[w_g-1]==0)  //the gene set becomes empty (i.e. gene g is the one and only gene in the gene set g)
		{
//			printf("\n Oh yes!!!");

			for (g=0; g < n_gene; g++) {
				if(w[g] > w_g) w[g] = w[g]-1;
			}
			
			for (s = w_g; s < n_gene; s++) {
				p[s-1] = p[s];
			}
			p[n_gene-1] = 0;
			
//			printf("p[0]=%d, p[1]=%d, p[2]=%d", p[0], p[1], p[2]);
			
			
			for (s = w_g; s < n_gene; s++) {
				for (i=0; i < n_sub; i++) {
					z[(s-1)*n_sub + i] = z[s*n_sub + i];
					n[(s-1)*n_sub + i] = n[s*n_sub + i];
					c[(s-1)*n_sub + i] = c[s*n_sub + i];
				}
			}

			for (i=0; i < n_sub; i++) {
				z[(n_gene-1)*n_sub + i] = 0;
				n[(n_gene-1)*n_sub + i] = 0;
				c[(n_gene-1)*n_sub + i] = 0;
			}
				
			*S = *S - 1;

			for (s = w_g; s < n_gene; s++) {
				K[s-1] = K[s];
				m[s-1] = m[s];
			}
			K[n_gene-1] = 0;
			m[n_gene-1] = 0;
		}	//if(p[w_g]=0)  //the gene set becomes empty (i.e. gene g is the one and only gene in the gene set g) 

	
		return;
	}

	double fn_lgamma(double x)
	{
		double x0,x2,xp,gl,gl0;
		int n,k;
		static double a[] = {
			8.333333333333333e-02,
			-2.777777777777778e-03,
			7.936507936507937e-04,
			-5.952380952380952e-04,
			8.417508417508418e-04,
			-1.917526917526918e-03,
			6.410256410256410e-03,
			-2.955065359477124e-02,
			1.796443723688307e-01,
			-1.39243221690590};
		
		x0 = x;
		if (x <= 0.0) return 1e308;
		else if ((x == 1.0) || (x == 2.0)) return 0.0;
		else if (x <= 7.0) {
			n = (int)(7-x);
			x0 = x+n;
		}
		
		x2 = 1.0/(x0*x0);
		xp = 2.0*M_PI;
		gl0 = a[9];
		
		for (k=8;k>=0;k--) {
			gl0 = gl0*x2 + a[k];
		}
		
		gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
		if (x <= 7.0) {
			for (k=1;k<=n;k++) {
				gl -= log(x0-1.0);
				x0 -= 1.0;
			}
		}
		return gl;
	}	


	//WHEN P==1, SAMPLE THETA_STAR
	void fn_theta_star_sampling_p1(double mu0, double tau20, double sig2, double y[], int n_sub, int c_is, int c[], int n, double *theta_star)
	{
		int j, i_sub;
		double var_tmp, mean_tmp, y_sum;
		
		var_tmp = 1.0/(1.0/tau20 + 1.0*n/sig2);
			
		y_sum= 0.0;
		for (i_sub=0; i_sub < n_sub; i_sub++) {
			if(c[i_sub]==c_is) y_sum = y_sum + y[i_sub];
		}
		mean_tmp = var_tmp*(mu0/tau20 + y_sum/sig2); 
			
		GetRNGstate();
		theta_star[c_is - 1] = rnorm(mean_tmp, sqrt(var_tmp));
		PutRNGstate();

		return;
	}
	
	
	//WHEN P==1, IMPUTE Z AND C FOR THE FIRST ITERATION.
	void fn_impute_z_c_iter1_p1(double M, double sig2, double mu0, double tau20, double pi_1, int n_sub, double y[], int *z, int *m, int *n, int *c, int *K)
	{
		int c_is, i_sub, k, j, i;
		double u, log_prob_max, tmp_prob;
		
		double *theta_star;
		double *prob_i;
		double *cumsum;
		
		//Allocate the memory
		theta_star =  new double[n_sub];  //MAX NUMBER OF CLUSTERS IS N_GENE
		prob_i =  new double[n_sub + 1];  //MAX NUMBER OF CLUSTERS IS N_GENE
		cumsum = new double[n_sub + 1]; 
		
		
		//INITIALIZATION
		for (i_sub=0; i_sub < n_sub; i_sub++) {
			z[i_sub] = 0;
			c[i_sub] = 0;
			n[i_sub] = 0;
		}
		
		
		//the first subject takes z_is=1 and makes the first cluster by itself
		z[0] = 1; *m = 1;
		*K = 1;
		n[0] = 1; c[0] = 1;
		
		c_is = c[0];
		fn_theta_star_sampling_p1(mu0, tau20, sig2, y, n_sub, c_is, c, n[c_is], & theta_star[0]);
		
		for(i_sub = 1; i_sub < n_sub ; i_sub++)
		{
			GetRNGstate();
			u = runif(0.0,1.0);
			PutRNGstate();
			
			if(u < pi_1) z[i_sub] = 1;
			
			if(z[i_sub]==1)
			{
				*m = *m + 1;
				
				for(k = 0; k < (*K); k++)
				{
					GetRNGstate();
					prob_i[k] = dnorm(y[i_sub], theta_star[k], sqrt(sig2), 1);
					PutRNGstate();
					
					prob_i[k] = log(1.0*n[k]/(M + 1.0*(*m) - 1.0)) + prob_i[k];
				}
				
				GetRNGstate();
				prob_i[*K] = dnorm(y[i_sub], mu0, sqrt(sig2 + tau20), 1);
				PutRNGstate();

				prob_i[*K] = log(M/(M + 1.0*(*m) - 1.0)) + prob_i[*K];
				
				
				log_prob_max=R_NegInf;
				for (k=0; k < (*K+1); k++) {
					if(log_prob_max < prob_i[k]) log_prob_max = prob_i[k];
				}
				
				// normalize the probabilities
				cumsum[0] = exp(prob_i[0] - log_prob_max);
				for(k=1; k < (*K+1); k++) cumsum[k] = cumsum[k-1] + exp(prob_i[k] - log_prob_max);
				for(k=0; k < (*K+1); k++) cumsum[k] = cumsum[k]/cumsum[*K];
				
				// Sample the cluster
				GetRNGstate();
				u = runif(0.0,1.0);
				PutRNGstate();
				
				for(k=0; k < (*K); k++) if(u <= cumsum[k]) break;
				c_is = k + 1;
				c[i_sub] = c_is;
				
				if(c_is == (*K+1))
				{
					//#START A NEW SUBJECT CLUSTER
					n[c_is - 1] = 1;
					fn_theta_star_sampling_p1(mu0, tau20, sig2, y, n_sub, c_is, c, n[c_is-1], & theta_star[0]);
					*K = *K + 1;
				}else{
					//#JOIN ONE OF EXISTING CLUSTERS
					n[c_is - 1] = n[c_is - 1] + 1;
					fn_theta_star_sampling_p1(mu0, tau20, sig2, y, n_sub, c_is, c, n[c_is-1], & theta_star[0]);
				} //if(c_is > (*K))
			} //if(z[i_sub]=1)
	
			//printf("\n i_sub=%d", i_sub);
			//printf("\n K_new = %d, m = %d", *K, *m);
			//for (i = 0; i < n_sub; i ++) {
			//	printf("\n i_sub=%d, z = %d, n= %d, c=%d", i, z[i], n[i], c[i]);
			//}
			
			
		} //for(i_sub = 1; i_sub < n_sub ; i_sub++)

//		printf("\n");
//		for(i_sub = 0; i_sub < n_sub ; i_sub++)
//		{
//			printf("c[%d]=%d,  ", i_sub, c[i_sub]);
//		}
		
		//Delete the memory
		delete[] theta_star; theta_star=  NULL;                
		delete[] prob_i; prob_i = NULL;  
		delete[] cumsum; cumsum = NULL;  
		
		return;
	}
	
	
	double fn_w_update_calc_prob (double pi_0, int m, double sig2, double tau20, double mu0, int z[], int n_sub, double tau21, double mu1, int n[], int K, int c[], double y[], double p, double Be, int G_p)
	{
		double prob, tmp_prob, y_s_sum, y_sum, tmp_var;
		int i_sub, i_k;
		
		prob = log(pi_0) + log(p/(Be + 1.0*G_p));
//		printf("\n Step1 prob=%f", prob);
		prob = prob - 1.0/2.0*m*log(2.0*M_PI*sig2) - 1.0/2.0*K*log(2*M_PI*tau20);
//		printf("\n Step2 prob=%f", prob);

		y_s_sum = 0.0;
		for (i_sub = 0; i_sub < n_sub; i_sub ++) {
			if(z[i_sub]==1) y_s_sum = y_s_sum + pow(y[i_sub], 2)/2.0/sig2;
		}

		prob = prob - 1.0*K*(pow(mu0, 2)/2.0/tau20) - y_s_sum;
//		printf("\n Step3 prob=%f", prob);

		prob = prob - 1.0*(n_sub - m)/2.0*log(2*M_PI*(tau21 + sig2));
//		printf("\n Step4 prob=%f", prob);
		
		if(m > 0) //#if any subject takes gene set s,
		{
			for(i_k = 0; i_k < K ; i_k ++) 
			{
				tmp_var = 1.0/(1.0/tau20 + n[i_k]/sig2);

				y_sum = 0.0;
				for (i_sub = 0; i_sub < n_sub; i_sub ++) {
					if(c[i_sub]==(i_k+1)) y_sum = y_sum + y[i_sub];
				}
				
//				printf("\n tmp_var=%f, y_sum=%f", tmp_var, y_sum);
				prob = prob + 1.0/2.0*log(2.0*M_PI*tmp_var) + 1.0/2.0*tmp_var*pow((mu0/tau20 + y_sum/sig2), 2); 
//				printf("\n Step5 i_k=%d, prob=%f", i_k, prob);
			}
		}
		
		if(m < n_sub)  //#if there is any subject not taking that gene set.
		{
			tmp_prob = 0.0;
			for (i_sub = 0; i_sub < n_sub; i_sub ++) {
				if(z[i_sub]==0) tmp_prob = tmp_prob + pow((mu1-y[i_sub]), 2)/2.0/(tau21 + sig2);
			}

			prob = prob - tmp_prob;
//			printf("\n Step6 prob=%f", prob);
		}
	
		return(prob);
	}

	void fn_remove_subject(int id_sub, int *z, int *m, int *c, int *k, int *n, int n_sub)
	{
		int cur_z, cur_c, i;
		
//		for (i=0; i < n_sub; i++) {
//			printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
//		}
//		printf("\n k=%d, m=%d", *k, *m);
		
		
		cur_z = z[id_sub];
		z[id_sub] = (-1);
		cur_c = c[id_sub];
		c[id_sub] = (-1);
		
//		printf("\n cur_z=%d", cur_z);
		if(cur_z == 1)
		{
			*m = (*m) - 1;
			n[cur_c-1] = n[cur_c-1] - 1;
			
			//IF A SUBJECT IS CLUSTERED BY ITSELF
			if (n[cur_c-1]==0) {
//				printf("\n Oh yes");
				for (i=0; i < n_sub; i++) {
					if (c[i] > cur_c) {
						c[i] = c[i] - 1;
					}
				}
				
				for (i=cur_c; i < n_sub; i++) {
					n[i-1] = n[i];
				}
				n[(n_sub-1)] = 0;
				*k = (*k) - 1;
			}
		}
		
//		for (i=0; i < n_sub; i++) {
//			printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
//		}
//		printf("\n k=%d, m=%d", *k, *m);
	
		return;
	}

	//update z and c for a given one gene set.
	void fn_z_c_one_gs_p1(double mu_0, double mu_1, double tau2_0, double tau2_1, double sig2, double pi_1, double M, int *z, int *m, int *c, int *K, int *n, int n_sub, double y[])
	{
		//save the current states
		double y_other, y_i;
		double *prob;
		double *cumsum;
		
		//Allocate the memory
		prob =  new double[(n_sub+2)];  //MAX NUMBER OF CLUSTERS IS N_SUB
		cumsum = new double[(n_sub+2)];  //MAX NUMBER OF CLUSTER IS N_SUB
		
		double log_prob_max, tmp_prob, tmp_1, tmp_2, tmp_var_1, tmp_var_2;
		double u;
		
		int i_sub, i, j, k;
		
		//		printf("p=%d, n_sub=%d", p, n_sub);
		
		
		for (i_sub=0; i_sub < n_sub; i_sub++)
		{
			//printf("\n i_sub=%d", i_sub);
			//for (i=0; i < n_sub; i++) {
			//printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
			//}
			//printf("\n k=%d, m=%d", *K, *m);
			
			fn_remove_subject(i_sub, z, m, c, K, n, n_sub);
			
			//for (i=0; i < n_sub; i++) {
			//printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
			//}
			//printf("\n k=%d, m=%d", *K, *m);
			
			y_i = y[i_sub];
			//printf("\n y_i[j]=%f", y_i[j]);
			
			//IF THERE IS AN EXISTING CLUSTER, THEN SUBJECT i MAY JOIN>
			if (*K > 0) {
				for (k=0; k < *K; k++) {
					
					y_other = 0.0;
					
					for (i=0; i < n_sub; i++) {
						if (c[i]==(k+1)) {
							y_other = y_other + y[i];
						}
					}
					
					tmp_prob = log(pi_1) + log(1.0*n[k]/(M + 1.0*(*m)));
					
					tmp_var_1 = 1.0/(1.0/tau2_0 + n[k]*1.0/sig2); 
					tmp_var_2 = 1.0/(1.0/tau2_0 + (n[k]*1.0 + 1.0)/sig2);
					tmp_1 = -1.0/2.0*log(tmp_var_1) - 1.0/2.0*tmp_var_1*pow((mu_0/tau2_0 + y_other/sig2), 2);
					tmp_2 = 1.0/2.0*log(tmp_var_2) + 1.0/2.0*tmp_var_2*pow((mu_0/tau2_0 + y_i/sig2 + y_other/sig2), 2);
						
					tmp_prob = tmp_prob - (1.0/2.0)*log(2.0*M_PI*sig2) - pow(y_i, 2)/2.0/sig2 + tmp_1 + tmp_2;

					prob[k] = tmp_prob;
				} //for (k=0; k < *K; k++)
			} //if (*K > 0) {
			
			//START A NEW CLUSTER OR BE AN INACTIVE SUBJECT
			tmp_1 = 0.0; tmp_2 = 0.0;
			tmp_1 = tmp_1 - (1.0/2.0)*log(2.0*M_PI*(tau2_0 + sig2)) - pow((mu_0-y_i), 2)/2.0/(tau2_0 + sig2); //START A NEW CLUSTER
			tmp_2 = tmp_2 - (1.0/2.0)*log(2.0*M_PI*(tau2_1 + sig2)) - pow((mu_1-y_i), 2)/2.0/(tau2_1 + sig2); //BE AN INACTIVE SUBJECT
			
			//START A NEW CLUSTER
			prob[*K] = log(pi_1) + log(M/(M+1.0*(*m))) + tmp_1;
			
			//BE AN INACTIVE SUBJECT
			prob[*K+1] = log(1-pi_1) + tmp_2;
			
			//printf("\n K= %d ", *K);
			//printf("\n ");
			//for (k=0; k < (*K+2); k++) {
			//printf("  prob[%d]=%f", k, prob[k]);
			//}
			//printf("\n");
			
			
			/////////////////////////////////////////////////////////////////////////////////////
			// normalize the probabilities
			/////////////////////////////////////////////////////////////////////////////////////
			log_prob_max=R_NegInf;
			for (k=0; k < (*K+2); k++) {
				if(log_prob_max < prob[k]) log_prob_max = prob[k];
			}
			
			cumsum[0] = exp(prob[0] - log_prob_max);
			for(k=1; k < (*K+2); k++) cumsum[k] = cumsum[k-1] + exp(prob[k] - log_prob_max);
			for(k=0; k < (*K+2); k++) cumsum[k] = cumsum[k]/cumsum[*K+1];
			
			/////////////////////////////////////////////////////////////////////////////////////
			// Sample the cluster
			/////////////////////////////////////////////////////////////////////////////////////
			GetRNGstate();
			u = runif(0.0,1.0);
			PutRNGstate();
			
			for(k=0; k < (*K+2); k++) if(u <= cumsum[k]) break;
			
			//printf("\n selected k= %d", k);
			
			if(k ==(*K+1))  //NOT TAKING A GENE SET S
			{
				z[i_sub] = 0;
				c[i_sub] = 0;
			}else{
				//TAKING A GENE SET S
				z[i_sub] = 1;
				*m = *m + 1;
				c[i_sub] = k+1;  //INDEX STARTS FROM 0 AND 0 STANDS FOR "INACTIVE" SO ADD 1 TO K
				
				if(k < *K)
				{
					//JOIN ONE EXISTING SUBJECT CLUSTER
					n[k] = n[k] + 1;
				}else{
					//START A NEW SUBJECT CLUSTER
					n[k] = 1;
					*K = *K + 1;
				}
			} //if(k =(*K+1))  #NOT TAKING A GENE SET S
			
			//printf("\n c_is=%d", k);
			//for (i=0; i < n_sub; i++) {
			//printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
			//}
			//printf("\n k=%d, m=%d", *K, *m);
			
		}  //for (i_sub=0; i_sub < n_sub; i_sub++)
		
		//		printf("\n Done!!!!");
		
		//Delete the memory
		delete[] prob; prob = NULL;  
		delete[] cumsum; cumsum =  NULL;                   
		
		return;
	}
	
	

	void fn_w_update(double *pi_0, double *pi_1, double *Be, double *M,
					 double *mu_0, double *tau2_0, double *mu_1, double *tau2_1, double *mu_2, double *tau2_2,
					 double *sig2,
					 int *w, int *p, int *G_p, int *S,
					 int *z, int *n, int *c, int *K, int *m,
					 double *y, int *nn_sub, int *nn_gene)
	{
		
//		double sig2 = *ssig2;
		int n_sub = *nn_sub;
		int n_gene = *nn_gene;
		
		int g, i_sub, i_gs, j, s, g_gs, S_old, k, i_g, i_iter; 
		int w_g, m_gs, K_gs;
		double tmp_prob, log_prob_max, u;
		double tau20_g, tau21_g, tau22_g, mu0_g, mu1_g, mu2_g, sig2_g;
		
		double *prob_gs;
		double *cumsum;
		double *y_g;
		int *z_gs;
		int *n_gs;
		int *c_gs;
		
		//Allocate the memory
		prob_gs =  new double[n_gene+2];  //MAX NUMBER OF CLUSTERS IS N_GENE
		cumsum = new double[n_gene+2]; //MAX NUMBER OF CLUSTER IS N_GENE
		y_g = new double[n_sub];
		z_gs = new int[n_sub];
		n_gs = new int[n_sub];
		c_gs = new int[n_sub];
		
		for (g = 0; g < n_gene; g++)
		{
//			printf("\n \n g=%d", g);
//			printf("\n");
//			for (i_g = 0; i_g < n_gene; i_g ++) {
//				printf("w[%d] = %d, ", i_g, w[i_g]);
//			}

//			printf("\n S=%d", *S);
//			for (i_gs = 0; i_gs < *S; i_gs ++)
//			{
//				printf("\n i_gs=%d", i_gs);
//				for (i_sub = 0; i_sub < n_sub; i_sub ++) {
//					printf("\n i_sub=%d, z = %d, n= %d, c=%d", i_sub, z[i_gs*n_sub + i_sub], n[i_gs*n_sub + i_sub], c[i_gs*n_sub + i_sub]);
//				}
//				printf("\n p=%d, m=%d, K=%d", p[i_gs], m[i_gs], K[i_gs]);
//			} //for (i_gs = 0; i_gs < *S; i_gs ++)

			//REMOVE ONE GENE G FROM GENE SET.
			w_g = w[g];
			w[g] = -1; 
			
//			printf("\n w_g=%d, p_g=%d", w_g, p[w_g-1]);
			if(w_g > 0)  //if a gene g is in a gene set.
			{
				fn_remove_gene(G_p, p, w, w_g, S, z, n, c, K, m, n_gene, n_sub);
			}
			
//			printf("\n");
//			for (i_g = 0; i_g < n_gene; i_g ++) {
//				printf("w[%d] = %d, ", i_g, w[i_g]);
//			}
//			printf("\n S=%d", *S);

//			for (i_gs = 0; i_gs < *S; i_gs ++)
//			{
//				printf("\n i_gs=%d", i_gs);
//				for (i_sub = 0; i_sub < n_sub; i_sub ++) {
//					printf("\n i_sub=%d, z = %d, n= %d, c=%d", i_sub, z[i_gs*n_sub + i_sub], n[i_gs*n_sub + i_sub], c[i_gs*n_sub + i_sub]);
//				}
//				printf("\n p=%d, m=%d, K=%d", p[i_gs], m[i_gs], K[i_gs]);
//			} //for (i_gs = 0; i_gs < *S; i_gs ++)
			
			
			//EXTRACT THINGS RELATED TO GENE G
			tau20_g = tau2_0[g]; tau21_g = tau2_1[g]; tau22_g = tau2_2[g];
			mu0_g = mu_0[g]; mu1_g = mu_1[g]; mu2_g = mu_2[g];
			sig2_g = sig2[g];
			
			for (i_sub = 0; i_sub < n_sub; i_sub ++) {
				y_g[i_sub] = y[g*n_sub + i_sub];
			}
			
//			for (i_sub = 0; i_sub < n_sub; i_sub ++) {
//				printf("\n y_g[%d] = %f", i_sub, y_g[i_sub]);
//			}
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//#JOINING ONE OF THE EXISTING GENE SETS
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(*S > 0)  //if there is any gene set (active + inactive)
			{
				for (i_gs = 0; i_gs < *S; i_gs ++)
				{
//					printf("\n i_gs=%d", i_gs);
					for (i_sub = 0; i_sub < n_sub; i_sub ++) {
						z_gs[i_sub] = z[i_gs*n_sub + i_sub];
						n_gs[i_sub] = n[i_gs*n_sub + i_sub];
						c_gs[i_sub] = c[i_gs*n_sub + i_sub];
						//						printf("\n i_sub=%d, z = %d, n= %d, c=%d", i_sub, z_gs[i_sub], n_gs[i_sub], c_gs[i_sub]);
					}
					//IF i_gs IS INACTIVE, m, z_gs, n_gs, K, c_gs ARE 0 (OR 0 VECTOR)
//					printf("\n m=%d, K=%d", m[i_gs], K[i_gs]);
					prob_gs[i_gs] = fn_w_update_calc_prob(*pi_0, m[i_gs], sig2_g, tau20_g, mu0_g, z_gs, n_sub, tau21_g, mu1_g, n_gs, K[i_gs], c_gs, y_g, (1.0*p[i_gs]), *Be, *G_p);
//					printf("\n i_gs=%d, prob=%f", i_gs, prob_gs[i_gs]);
				} //for (i_gs = 0; i_gs < S_p; i_gs ++)
			} //if(S_p > 0)  //if there is any gene set (active + inactive)
			
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//#STARTING A NEW GENE SET BY ITSELF
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//#IMPUTE Z AND C FOR STARTING A NEW GENE SET
			fn_impute_z_c_iter1_p1(*M, sig2_g, mu0_g, tau20_g, *pi_1, n_sub, y_g, &z_gs[0], &m_gs, &n_gs[0], &c_gs[0], &K_gs);

			//printf("\n K_new = %d, m = %d", K_gs, m_gs);
			//for (i_sub = 0; i_sub < n_sub; i_sub ++) {
			//	printf("\n i_sub=%d, z = %d, n= %d, c=%d", i_sub, z_gs[i_sub], n_gs[i_sub], c_gs[i_sub]);
			//}
			
			//CALL 100 TIMES FOR BETTER MIXING	
			for (i_iter=0; i_iter < 30 ; i_iter ++) {
				fn_z_c_one_gs_p1(mu0_g, mu1_g, tau20_g, tau21_g, sig2_g, *pi_1, *M, &z_gs[0], &m_gs, &c_gs[0], &K_gs, &n_gs[0], n_sub, y_g);
			//	printf("\n K_new = %d", K_gs);
			}
			
			tmp_prob = fn_w_update_calc_prob(*pi_0, m_gs, sig2_g, tau20_g, mu0_g, z_gs, n_sub, tau21_g, mu1_g, n_gs, K_gs, c_gs, y_g, *Be, *Be, *G_p);
			
//			printf("\n tmp_prob=%f", tmp_prob);
			if(K_gs > 0)
			{
				tmp_prob = tmp_prob + 1.0*m_gs*log(*pi_1) + 1.0*(n_sub - m_gs)*log(1 - *pi_1) + 1.0*K_gs*log(*M); 
				
				for (k=0; k < K_gs; k++) {
					tmp_prob = tmp_prob + fn_lgamma(1.0*n_gs[k]);
				}
				
				for (j=0; j < m_gs; j++) {
					tmp_prob = tmp_prob - log(*M + j);
				}
			}else {
				tmp_prob = tmp_prob + 1.0*n_sub*log(1 - *pi_1);
			}
			
			prob_gs[*S] = tmp_prob;
//			printf("\n start a new cluster, prob=%f", prob_gs[*S + 1]);
			
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//#remain as an inactive gene
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			tmp_prob = log(1.0 - *pi_0) - n_sub/2.0*log(2.0*M_PI*(tau22_g + sig2_g)); 
			for (i_sub = 0; i_sub < n_sub; i_sub++) {
				tmp_prob = tmp_prob - pow((mu2_g - y_g[i_sub]), 2)/2.0/(tau22_g + sig2_g);
			}
			prob_gs[*S + 1] = tmp_prob;
//			printf("\n remain inactive, prob=%f", prob_gs[*S_p + 2]);
			
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//sample one gs
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//			printf("\n");
//			for (s=0; s < (*S+2); s++) {
//				printf("prob[%d]=%f  ", s, prob_gs[s]);
//			}
			
			
			log_prob_max=R_NegInf;
			for (s=0; s < (*S+2); s++) {
				if(log_prob_max < prob_gs[s]) log_prob_max = prob_gs[s];
			}
			
			// normalize the probabilities
			cumsum[0] = exp(prob_gs[0] - log_prob_max);
			for (s = 1; s < (*S+2); s++) cumsum[s] = cumsum[s-1] + exp(prob_gs[s] - log_prob_max);
			for (s = 0; s < (*S+2); s++) cumsum[s] = cumsum[s]/cumsum[*S+1];
			
			// Sample the cluster
			GetRNGstate();
			u = runif(0.0,1.0);
			PutRNGstate();
			
			for(s = 0; s < (*S+2); s++) if(u <= cumsum[s]) break;
			g_gs = s + 1;  //g_gs=0 implies inactive, so gene cluster is indexed from 1 to S^\prime
			
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//#UPDATE THE CONFIGURATION ACCORDING TO THE SAMPLE OF w_g
			//#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			S_old = *S;
			
			if(g_gs <= S_old) //#joining an existing gene set
			{
				w[g] = g_gs;
				p[g_gs-1] = p[g_gs-1] + 1;	
				*G_p = *G_p + 1;
				//#gam.star and theta will be updated later
			}
			
			if(g_gs == (S_old + 1))  //#start a new gene set.
			{
				*G_p = *G_p + 1;
				*S = *S + 1; 
				
				w[g] = *S; 
				
				for (i_sub = 0; i_sub < n_sub; i_sub++) {
					z[(g_gs-1)*n_sub + i_sub] = z_gs[i_sub];
					c[(g_gs-1)*n_sub + i_sub] = c_gs[i_sub];
					n[(g_gs-1)*n_sub + i_sub] = n_gs[i_sub];
				}
				
				K[g_gs-1] = K_gs;
				m[g_gs-1] = m_gs;
				p[g_gs - 1] = 1;
			} //if(g_gs = (S_old + 1))  //#start a new gene set.
			
			if(g_gs == (S_old + 2)) { w[g]= 0; } // #an inactive gene
			
//			printf("\n g=%d, g_gs=%d, w[g]=%d", g, g_gs, w[g]);
			
			
//			printf("\n");
//			for (i_g = 0; i_g < n_gene; i_g ++) {
//				printf("w[%d] = %d, ", i_g, w[i_g]);
//			}
			
//			for (i_gs = 0; i_gs < *S; i_gs ++)
//			{
//				printf("\n i_gs=%d", i_gs);
//				for (i_sub = 0; i_sub < n_sub; i_sub ++) {
//					printf("\n i_sub=%d, z = %d, n= %d, c=%d", i_sub, z[i_gs*n_sub + i_sub], n[i_gs*n_sub + i_sub], c[i_gs*n_sub + i_sub]);
//				}
//				printf("\n p=%d, m=%d, K=%d", p[i_gs], m[i_gs], K[i_gs]);
//			} //for (i_gs = 0; i_gs < S; i_gs ++)
			
			
		} //for (g = 0; g < n_gene; g++) {
		
		//Delete the memory
		delete[] prob_gs; prob_gs=  NULL;                
		delete[] cumsum; cumsum = NULL;  
		delete[] y_g; y_g =  NULL;                
		delete[] z_gs; z_gs =  NULL;                
		delete[] n_gs; n_gs =  NULL;                
		delete[] c_gs; c_gs =  NULL;                

		return; 
	
	}
	
}




