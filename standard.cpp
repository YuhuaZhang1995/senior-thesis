#include "header.h"
#include <cmath>
#include <iostream>

using namespace std;

void LD_mat(double* sigma, double* X, int snpnum, int samplesize){
//    double* sigma=new double[snpnum*snpnum];
    int i;
    int j;
    int m;
    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            for(m=0;m<samplesize;m++){
                sigma[i*snpnum+j]+=X[i*samplesize+m]*X[j*samplesize+m];
            }
        }
    }
}

void cal_corr(double* sigma,int snpnum, double* sigma_k, double* CorrelationMatrix, double sigma_e, int samplesize){
    double* temp=new double[snpnum*snpnum];
    int i;
    int j;
    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
//            temp[i*snpnum+j]=0;
            temp[i*snpnum+j]=sigma[i*snpnum+j];
        }
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            temp[j*snpnum+i]=temp[j*snpnum+i]*sigma_k[i];
//            cout<<"temp"<<temp[j*snpnum+i]<<"\t";
        }
//        cout<<"\n";
    }
    int m;
    double *tmp=new double[snpnum*snpnum];

    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            tmp[i*snpnum+j]=0;
        }
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            for(m=0;m<snpnum;m++){
                tmp[i*snpnum+j]=tmp[i*snpnum+j]+temp[i*snpnum+m]*sigma[m*snpnum+j];
            }
//            cout<<"tmp"<<tmp[i*snpnum+j]<<"\t";
        }
//        cout<<"\n";
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            tmp[i*snpnum+j]=tmp[i*snpnum+j]*samplesize/sigma_e+sigma[i*snpnum+j];
//            cout<<"tmp:"<<tmp[i*snpnum+j]<<"\t";
        }
//        cout<<"\n";
    }

    for(i=1;i<snpnum;i++){
        for(j=0;j<i;j++){
            CorrelationMatrix[(i*(i-1))/2+j]=tmp[i*snpnum+j]/sqrt(tmp[i*snpnum+i]*tmp[j*snpnum+j]);
//            cout<<CorrelationMatrix[(i*(i-1))/2+j]<<"\t";
        }
//        cout<<"\n";
    }

    delete temp;
    delete tmp;
}
