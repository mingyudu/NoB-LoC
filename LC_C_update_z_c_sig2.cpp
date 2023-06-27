#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>

extern "C"
{
//#include "LC_C_update_w.h"

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
	void fn_z_c_one_gs(double *mu_0, double *mu_1, double *tau2_0, double *tau2_1, double *sig2, double *ppi_1, double *MM, int *z, int *m, int *c, int *K, int *n, int *nn_sub, int*pp, double *y)
	{
//		double sig2 = *ssig2;
		double pi_1 = *ppi_1;
		double M = *MM;
		int n_sub = *nn_sub;
		int p = *pp;
		
		//save the current states
		double *y_other;
		double *y_i;
		double *prob;
		double *cumsum;
		
		//Allocate the memory
		y_other =  new double[p]; 
		y_i =  new double[p]; 
		prob =  new double[(n_sub+2)];  //MAX NUMBER OF CLUSTERS IS N_SUB
		cumsum = new double[(n_sub+2)];  //MAX NUMBER OF CLUSTER IS N_SUB
		
		double log_prob_max, tmp_prob, tmp_1, tmp_2, tmp_var_1, tmp_var_2;
		double u;
		
		int i_sub, i, j, k;
		
//		printf("p=%d, n_sub=%d", p, n_sub);

		
		for (i_sub=0; i_sub < n_sub; i_sub++)
		{
//			printf("\n i_sub=%d", i_sub);
//			for (i=0; i < n_sub; i++) {
//				printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
//			}
//			printf("\n k=%d, m=%d", *K, *m);
			
			fn_remove_subject(i_sub, z, m, c, K, n, n_sub);
			
//			for (i=0; i < n_sub; i++) {
//				printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
//			}
//			printf("\n k=%d, m=%d", *K, *m);
			
			
			for (j=0; j < p; j++) {
				y_i[j] = y[j*n_sub + i_sub];
//				printf("\n y_i[j]=%f", y_i[j]);
			}
			
			//IF THERE IS AN EXISTING CLUSTER, THEN SUBJECT i MAY JOIN>
			if (*K > 0) {
				for (k=0; k < *K; k++) {
					
					for (j=0; j < p ; j++) {
						y_other[j] = 0.0;
					}
					
					for (i=0; i < n_sub; i++) {
						if (c[i]==(k+1)) {
							for (j=0; j < p ; j++) {
								y_other[j] = y_other[j] + y[j*n_sub + i];
							}
						}
					}
					
					tmp_prob = log(pi_1) + log(1.0*n[k]/(M + 1.0*(*m)));
					
					for (j=0; j < p; j++) {
						tmp_var_1 = 1.0/(1.0/tau2_0[j] + n[k]*1.0/sig2[j]); 
						tmp_var_2 = 1.0/(1.0/tau2_0[j] + (n[k]*1.0 + 1.0)/sig2[j]);
						tmp_1 = -1.0/2.0*log(tmp_var_1) - 1.0/2.0*tmp_var_1*pow((mu_0[j]/tau2_0[j] + y_other[j]/sig2[j]), 2);
						tmp_2 = 1.0/2.0*log(tmp_var_2) + 1.0/2.0*tmp_var_2*pow((mu_0[j]/tau2_0[j] + y_i[j]/sig2[j] + y_other[j]/sig2[j]), 2);
						
						tmp_prob = tmp_prob - (1.0/2.0)*log(2.0*M_PI*sig2[j]) - pow(y_i[j], 2)/2.0/sig2[j] + tmp_1 + tmp_2;
					}
					prob[k] = tmp_prob;
				} //for (k=0; k < *K; k++)
			} //if (*K > 0) {
			
			//START A NEW CLUSTER OR BE AN INACTIVE SUBJECT
			tmp_1 = 0.0; tmp_2 = 0.0;
			for (j=0; j < p; j++) {
				tmp_1 = tmp_1 - (1.0/2.0)*log(2.0*M_PI*(tau2_0[j] + sig2[j])) - pow((mu_0[j]-y_i[j]), 2)/2.0/(tau2_0[j] + sig2[j]); //START A NEW CLUSTER
				tmp_2 = tmp_2 - (1.0/2.0)*log(2.0*M_PI*(tau2_1[j] + sig2[j])) - pow((mu_1[j]-y_i[j]), 2)/2.0/(tau2_1[j] + sig2[j]); //BE AN INACTIVE SUBJECT
			}
			
			//START A NEW CLUSTER
			prob[*K] = log(pi_1) + log(M/(M+1.0*(*m))) + tmp_1;
			
			//BE AN INACTIVE SUBJECT
			prob[*K+1] = log(1-pi_1) + tmp_2;
			
//			printf("\n K= %d ", *K);
//			printf("\n ");
//			for (k=0; k < (*K+2); k++) {
//				printf("  prob[%d]=%f", k, prob[k]);
//			}
//			printf("\n");


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
			
//			printf("\n selected k= %d", k);
			
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
			
//			printf("\n c_is=%d", k);
//			for (i=0; i < n_sub; i++) {
//				printf("\n i= %d, z[i] = %d, c[i]=%d, n[i]=%d", i, z[i], c[i], n[i]);
//			}
//			printf("\n k=%d, m=%d", *K, *m);
			
		}  //for (i_sub=0; i_sub < n_sub; i_sub++)
		
//		printf("\n Done!!!!");
		
		//Delete the memory
		delete[] y_i; y_i =  NULL;                
		delete[] y_other; y_other =  NULL;                
		delete[] prob; prob = NULL;  
		delete[] cumsum; cumsum =  NULL;                   
		
		return;
	}
	
}




