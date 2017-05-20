#include "header.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//#include <vector>
#include <iostream>
//#include <array>

using namespace std;

void poisson_dist(double* gamma,int* anno, int k, double* pi,int s_anno){
    int i;
    double sum_anno=0;
    double temp1=1;
    double temp2=1;
    double sum=0;
//    double pk;

    for(i=0;i<s_anno;i++){
        sum_anno=sum_anno+gamma[i]*anno[i];
    }

    if(sum_anno<0){
        sum_anno=0.01;
    }

    for(i=0;i<=k;i++){
        if(i==0){
            pi[0]=1/(1*exp(sum_anno));
            sum=sum+pi[i];
        }
        else{
            for(i=1;i<=k;i++){
                temp1=sum_anno*temp1;
                temp2=temp2*i;
                pi[i]=temp1/(temp2*exp(sum_anno));
                sum=sum+pi[i];
            }
        }
    }
    sum_anno=0;
    temp1=1;
    temp2=1;
    //To make sure the summation of the probability id 1
    for(i=0;i<=k;i++){
        pi[i]=pi[i]/sum;
    }
    sum=0;
}

//to store the matrix information into an array, each k+1 equals to one snp's information
 void Calc_Priors(int k, int snpnum,/* double* anno,*/ double* gamma,int* indicator/*n*j*/, double* poi_prob/*n*(k+1)*/,int s_anno){
//    cout<<s_anno;
    int* anno_=new int[s_anno];
    double* temp=new double[k+1];
    int i;
    int j;
    for(i=0;i<snpnum;i++){
        for(j=0;j<s_anno;j++){
            anno_[j]=indicator[i*s_anno+j];
 //           cout<<"anno"<<anno_[j]<<"\t";
        }
 //       cout<<"\n";
        poisson_dist(gamma,anno_,k,temp,s_anno);
        for(j=0;j<=k;j++){
            poi_prob[i*(k+1)+j]=temp[j];
 //           cout<<"temp"<<temp[j]<<"\t";
        }
 //       cout<<"\n";
/*        for(j=0;j<s_anno;j++){
            delete anno_[j];
        }
        for(j=0;j<=k;j++){
            delete temp[j];
        }*/
    }
}
