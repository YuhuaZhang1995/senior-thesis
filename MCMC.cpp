#include "header.h"
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;

//sample data according to the poisson distribution
void MCMC_smaple(double* mar_prob, double* post_prob, double* z_score, int samplesize, double *sigma, int k, int snpnum, double* gamma, int*indicator, int s_anno){
    //return with the matrix to sample from
    double* poi_prob = new double[snpnum*(k+1)];
    Calc_Priors(k, snpnum,gamma,indicator,poi_prob, s_anno);


    int iter = 10;
    int* snp_order=new int[iter*snpnum];

    int i;
    int j;
    int m;
    double* temp_prob=new double[snpnum*(k+1)];
    for(i=0;i<snpnum;i++){
        temp_prob[i*(k+1)]=poi_prob[i*(k+1)];
        for(j=1;j<=k;j++){
            temp_prob[i*(k+1)+j]=temp_prob[i*(k+1)+j-1]+poi_prob[i*(k+1)+j];
//            cout<<"temp:"<<temp_prob[i*(k+1)+j]<<"\t";
        }
//        cout<<"\n";
    }

    double random_seed;
    srand(time(0));
    for(i=0;i<iter;i++){
        for(j=0;j<snpnum;j++){
            //generate random value
            random_seed = ((double)rand()/RAND_MAX);
            for(m=0;m<=k;m++){
                if(random_seed<=temp_prob[j*(k+1)+m]){
                    snp_order[i*snpnum+j]=m;
                    break;
                }
 //           cout<<"\n";
            }
//            cout<<"s_o"<<snp_order[i*snpnum+j]<<"\t";
//            cout<<"r:"<<random_seed<<"\t";
        }
//        cout<<"\n";
    }
/*    for(i=0;i<10;i++){
        for(j=0;j<snpnum;j++){
            cout<<snp_order[i*(k+1)+j]<<"\t";
        }
        cout<<"\n";
    }
*/

    double intergration=0;
    int* temp=new int[snpnum];
    double* sigma_k=new double[snpnum];
    double* CorrelationMatrix=new double[snpnum*(snpnum-1)/2];
    double sigma_e=1;
//    double temp_num;
    int index;

    for(i=0;i<iter;i++){

        for(j=0;j<snpnum;j++){
            temp[j]=snp_order[i*snpnum+j];
            sigma_k[j]=temp[j];//may have some problems here
        }

        cal_corr(sigma, snpnum, sigma_k, CorrelationMatrix, sigma_e, samplesize);
        intergration=pmvnorm(snpnum, z_score, CorrelationMatrix);
 //       cout<<"i:"<<intergration<<"\n";

/*        temp_num=1;
        for(j=0;j<snpnum;j++){
            index=temp[j];
            temp_num=temp_num*poi_prob[j*(k+1)+index];
        }
        temp_num=temp_num*intergration;
*/
        for(j=0;j<snpnum;j++){
            index=temp[j];
            post_prob[j*(k+1)+index] +=intergration;
        }
//        cout<<"p:"<<post_prob[i*(k+1)+j]<<"\n";
    }
    double sum=0;
    for(i=0;i<snpnum;i++){
        for(j=0;j<=k;j++){
            sum+=post_prob[i*(k+1)+j];
        }
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<=k;j++){
            post_prob[i*(k+1)+j]=post_prob[i*(k+1)+j]/sum;
            cout<<"pp:"<<post_prob[i*(k+1)+j]<<"\t";
        }
        cout<<"\n";
    }
    for(i=0;i<snpnum;i++){
        sum=0;
        for(j=0;j<=k;j++){
            sum+=post_prob[i*(k+1)+j];
        }
//        cout <<"sum:"<<sum<<"\n";
        mar_prob[i]=(sum-post_prob[i*(k+1)])/sum;
//        cout <<"mar:"<<mar_prob[i]<<"\n";
    }

    delete snp_order;
    delete temp;
    delete poi_prob;
}
