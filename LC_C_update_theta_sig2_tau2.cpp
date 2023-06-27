#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>

extern "C"
{
	//THETAS FOR INACTIVE GENES AND ALL PATIENTS
	void fn_theta_sub2(double y_g[], double sig2, double mu, double a_2, double b_2, double *tau2, int N_sub, double *theta, int gene_id)
	{
		int i_sub;
		double v, m;

		double sq_dev = 0.0;
		double a_tmp, b_tmp;

		for (i_sub=0; i_sub < N_sub; i_sub++) {
			v = 1.0/(1.0/tau2[gene_id] + 1.0/sig2);
			m = v*(mu/tau2[gene_id] + y_g[i_sub]/sig2);
				
			GetRNGstate();
			theta[gene_id*N_sub + i_sub] = rnorm(m, sqrt(v));
			PutRNGstate();

			sq_dev = sq_dev + pow((theta[gene_id*N_sub + i_sub] - mu), 2);
		}

		//UPDATE TAU2_G
		a_tmp = a_2 + 1.0*N_sub/2.0;
		b_tmp = b_2 + sq_dev/2.0;
		GetRNGstate();
		tau2[gene_id] = 1.0/rgamma(a_tmp, 1.0/b_tmp);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
		PutRNGstate();
		
		
		return;
	}

	//THETAS FOR ACTIVE GENES AND ALL PATIENTS
	void fn_theta_sub1(double y_g[], double sig2, double mu0, double a_0, double b_0, double *tau20, double mu1, double a_1, double b_1, double *tau21, int n[], int K, int c[], double *theta, int N_sub, int gene_id)
	{
		int i_c, i_sub;
		double theta_star, v, m, y_sum;
		
		double sq_dev_0, sq_dev_1;
		double a_tmp, b_tmp;
		int count_p;
		
		//////////////////////////////////////////////////////////////////
		//ACTIVE PATIENTS
		//////////////////////////////////////////////////////////////////
		sq_dev_0 = 0.0;
		for (i_c=1; i_c <= K; i_c++) {
			
			y_sum = 0.0;
			for (i_sub=0; i_sub < N_sub; i_sub++) {
				if(c[i_sub]==i_c) {y_sum = y_sum + y_g[i_sub];}
			}

			v= 1.0/(1.0/tau20[gene_id] + 1.0*n[i_c-1]/sig2);
			m=v*(mu0/tau20[gene_id] + y_sum/sig2);
			
			GetRNGstate();
			theta_star = rnorm(m, sqrt(v));
			PutRNGstate();
			
			for (i_sub=0; i_sub < N_sub; i_sub++) {
				if(c[i_sub]==i_c) {theta[gene_id*N_sub + i_sub] = theta_star;}
			}
			
			sq_dev_0 = sq_dev_0 + pow((theta_star-mu0), 2);
			
		} //for (i_c=0; i_c < K; i_c++) {

		
		//UPDATE TAU0_G
		a_tmp = a_0 + 1.0*K/2.0;
		b_tmp = b_0 + sq_dev_0/2.0;
		GetRNGstate();
		tau20[gene_id] = 1.0/rgamma(a_tmp, 1.0/b_tmp);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
		PutRNGstate();
		
		//////////////////////////////////////////////////////////////////
		//INACTIVE PATIENTS
		//////////////////////////////////////////////////////////////////
		count_p = 0; sq_dev_0 = 0.0;
		for (i_sub=0; i_sub < N_sub; i_sub++) {
			if(c[i_sub]==0){
				v= 1.0/(1.0/tau21[gene_id] + 1/sig2);
				m=v*(mu1/tau21[gene_id] + y_g[i_sub]/sig2);
				
				GetRNGstate();
				theta[gene_id*N_sub + i_sub] = rnorm(m, sqrt(v));
				PutRNGstate();
				
				count_p = 1 + count_p;
				sq_dev_1 = sq_dev_1 + pow((theta[gene_id*N_sub + i_sub] - mu1), 2);
			}
		}

		//UPDATE TAU1_G
		a_tmp = a_1 + 1.0*count_p/2.0;
		b_tmp = b_1 + sq_dev_1/2.0;
		GetRNGstate();
		tau21[gene_id] = 1.0/rgamma(a_tmp, 1.0/b_tmp);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
		PutRNGstate();
		
		return;
	}

	void fn_sig2_sub(double *a_g, double *b_g, int N_sub, int N_gene, double sq_dev[], double *sig2)
	{
		int i_g;
		double a_tmp, b_tmp;
		
		for (i_g=0; i_g < N_gene; i_g++) {
			a_tmp = a_g[i_g] + 1.0*N_sub/2.0;
			b_tmp = b_g[i_g] + sq_dev[i_g]/2.0;
			GetRNGstate();
			sig2[i_g] = 1.0/rgamma(a_tmp, 1.0/b_tmp);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
			PutRNGstate();
		}
		return;
	}				 
	
	void fn_theta_sig2_update(int *NN_gene, int *NN_sub, 
						 double *mu0, double  *tau20, double *mu1, double *tau21, double *mu2, double *tau22,
						 double *a_g, double *b_g,
						 double *a0_g, double *b0_g, double *a1_g, double *b1_g, double *a2_g, double *b2_g,     
						 double *y, double *theta, double *sig2,
						 int *w, int *n, int *c, int *K)  
	{
		int N_gene = *NN_gene;
		int N_sub = *NN_sub;

		int i_g, w_g, i_sub;
		
		double *y_g;
		int *n_g;
		int *c_g;
		double *sq_dev;
		
		//Allocate the memory
		y_g =  new double[N_sub];  
		n_g =  new int[N_sub];  
		c_g = new int[N_sub]; 
		sq_dev = new double[N_gene]; 
		
		for (i_g=0; i_g < N_gene; i_g++) {
			if (w[i_g]==0) {
				//INACTIVE GENES
				for(i_sub=0; i_sub < N_sub; i_sub ++)
				{
					y_g[i_sub] = y[i_g*N_sub + i_sub];
				}
				fn_theta_sub2(y_g, sig2[i_g], mu2[i_g], a2_g[i_g], b2_g[i_g], tau22, N_sub, theta, i_g);
				
				//UPDATE TAU0_G AND TAU1_G
				//SAMPLE FROM ITS PRIOR
				GetRNGstate();
				tau20[i_g] = 1.0/rgamma(a0_g[i_g], 1.0/b0_g[i_g]);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
				PutRNGstate();
				
				GetRNGstate();
				tau21[i_g] = 1.0/rgamma(a1_g[i_g], 1.0/b1_g[i_g]);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
				PutRNGstate();
			}else {
				//ACTIVE GENES
				w_g=w[i_g];
				
				for(i_sub=0; i_sub < N_sub; i_sub ++)
				{
					y_g[i_sub] = y[i_g*N_sub + i_sub];
					c_g[i_sub] = c[(w_g - 1)*N_sub + i_sub];
					n_g[i_sub] = n[(w_g - 1)*N_sub + i_sub];
				}
				fn_theta_sub1(y_g, sig2[i_g], mu0[i_g], a0_g[i_g], b0_g[i_g], tau20, mu1[i_g], a1_g[i_g], b1_g[i_g], tau21, n_g, K[w_g-1], c_g, theta, N_sub, i_g);
				
				//UPDATE TAU2_G
				//SAMPLE FROM ITS PRIOR
				GetRNGstate();
				tau22[i_g] = 1.0/rgamma(a2_g[i_g], 1.0/b2_g[i_g]);  //INVERSE-GAMMA (A, B) WITH MEAN=B/(A-1)
				PutRNGstate();
			}//if (w[i_g]==0) {

		} //for (i_g=0; i_g < N_gene; i_g++)

	
		//siq2 SAMPLING	
		for (i_g=0; i_g < N_gene; i_g++) {
			sq_dev[i_g] = 0.0;
			for (i_sub = 0; i_sub < N_sub; i_sub++) {
				sq_dev[i_g] = sq_dev[i_g] + pow(y[i_g*N_sub + i_sub] - theta[i_g*N_sub + i_sub],2);
			}
		}
		
		fn_sig2_sub(a_g, b_g, N_sub, N_gene, sq_dev, sig2);
		
		//Delete the memory
		delete[] y_g; y_g=  NULL;                
		delete[] n_g; n_g = NULL;  
		delete[] c_g; c_g=  NULL;                
	
		return;
	}
	
	
}




