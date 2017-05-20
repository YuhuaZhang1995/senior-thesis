#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

extern"C" {

extern void mvndst_(int* n, double* lower, double* upper, int* infin, double* correl, int* maxpts, double* abseps, double* releps, double* error, double* value, int* inform);
}

double mvnnorm(int* n, double* lower, double* upper, int* infin, double* correl, int* maxpts, double* abseps, double* releps, double* error, double* value, int* inform);
double pmvnorm(int n, double* bound, double* correlationMatrix);

void cal_corr(double* sigma,int snpnum, double* sigma_k, double* CorrelationMatrix, double sigma_e, int samplesize);

void poisson_dist(double* gamma, double* anno, int k, double* pi);
void Calc_Priors(int k, int snpnum, double* gamma,int* indicator/*n*(k+1)*/, double* poi_prob/*n*(k+1)*/, int s_anno);

//void Estep_build(int snpnum, int k, int* snp_order);
//void Estep_search(int s_anno, double* mar_prob, int samplesize,int k, int snpnum, double* z_score, double* CorrelationMatrix, double* gamma,int* indicator, double* post_prob);

void MCMC_smaple(double* mar_prob, double* post_prob, double* z_score, int samplesize, double *sigma, int k, int snpnum, double* gamma, int*indicator, int s_anno);

void Mstep_prepare(int anno_size, double* post_prob,int snpnum, int k, int* anno_, double* index/*anno_size*(snpnum+1)*/);
//int Mstep_f (const gsl_vector* x, void* params, gsl_vector* f);
void Mstep_solve(int snpnum, int k, int anno_size, int* anno_, double* post_prob, double* new_gamma);

#endif // HEADER_H_INCLUDED
