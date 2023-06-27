#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include <math.h>

extern "C"
{
    
void Calc_Pair_Prob(int *nn, int *TT, int *Ss_sample, int *Prob)
{
	int		n = *nn; 			// n observations
	int		T = *TT;			// n Gibbs draws without burn-in
		
	int i, j, k;		// loop
	int tmp_i, tmp_j;

	int *s_cur;                            //combination at the current state
	s_cur = new int[n];                       //cluster size at the current state
		
	for (k=0; k < T; k++)
	{
		for (i=0; i < n; i++) s_cur[i] = Ss_sample[k*n+i];
			
		for (i=0; i < (n-1); i++)
		{
			tmp_i = s_cur[i];
			//printf("\n i=%d, s[i]=%d", i, tmp_i);
			
			if (tmp_i > 0) {
				for (j=(i+1); j < n; j++) {
					tmp_j = s_cur[j];
					
					if (tmp_i==tmp_j) {
						Prob[n*i+j] = Prob[n*i+j] + 1;
					}
					//printf("\n i=%d, s[i]=%d, j=%d, s[j]=%d, Prob[n*i+j]=%d", i, tmp_i, j, tmp_j, Prob[n*i+j]);
					
				} //for (j=(i+1); j < n; j++)
			} //if (tmp_i > 0)
		} //for (i=0; i < (n-1); i++)
	} //for (k=0; k < T; k++)		
			
	//Delete the memory
	delete[] s_cur; s_cur =  NULL;                   
	return;
}	

	//FIND THE LEAST SQUARES CLUSTERING (DAVID DAHL (2006))
	void Calc_Squared_Distance(int *nn, int *TT, int *Ss_sample, double *Prob, double *sq_dis)
	{
		int		n = *nn; 			// n observations
		int		T = *TT;			// n Gibbs draws without burn-in
		
		int i, j, k;		// loop
		int tmp_i, tmp_j;
		double sq_dis_k = 0.0;
		
		int *s_cur;                            //combination at the current state
		s_cur = new int[n];                    //cluster size at the current state
		
		for (k=0; k < T; k++)
		{
			for (i=0; i < n; i++) s_cur[i] = Ss_sample[k*n+i];
			
			sq_dis_k = 0.0;
			for (i=0; i < (n-1); i++)
			{
				tmp_i = s_cur[i];
				//printf("\n i=%d, s[i]=%d", i, tmp_i);
				
				if (tmp_i > 0) {
					for (j=(i+1); j < n; j++) {
						tmp_j = s_cur[j];
						
						if (tmp_i==tmp_j) {
							sq_dis_k = sq_dis_k + pow((1 - Prob[n*i+j]), 2);
						}else {
							sq_dis_k = sq_dis_k + pow((0 - Prob[n*i+j]), 2);
						} //if (tmp_i==tmp_j) {
						//printf("\n i=%d, s[i]=%d, j=%d, s[j]=%d, Prob[n*i+j]=%d", i, tmp_i, j, tmp_j, Prob[n*i+j]);
					} //for (j=(i+1); j < n; j++)
				}else {
					for (j=(i+1); j < n; j++) {
						sq_dis_k = sq_dis_k + pow((0 - Prob[n*i+j]), 2);
					}
				}//if (tmp_i > 0)
			} //for (i=0; i < (n-1); i++)
				
			sq_dis[k] = sq_dis_k;
		} //for (k=0; k < T; k++)		
		
		//Delete the memory
		delete[] s_cur; s_cur =  NULL;                   
		return;
	}	
	
}
