#include "header.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;
//estimate the maximal value get from Estep
//currently, decide to use Newton method to evaluate the solver of the equation
//remember to install the GSL library

void Mstep_prepare(int anno_size, double* post_prob,int snpnum, int k, int* anno_, double* index/*anno_size*(snpnum+1)*/){
//    double* index=new double[anno_size*(snpnum+1)]; //store the parameters of the equations

    int i;
    int j;
    int m;
//    double temp;
/*    for(i=0;i<anno_size;i++){
        for(j=0;j<=snpnum;j++){
            index[i*(snpnum+1)+j]=0;
        }
    }*/

    for(i=0;i<anno_size;i++){
        for(j=0;j<snpnum;j++){
  //          temp=0;
            for(m=0;m<=k;m++){
                index[i*(snpnum+1)+j]+=post_prob[j*(k+1)+m]*m;
//               temp=temp-post_prob[j*(k+1)+m];
//               cout<<post_prob[j*(k+1)+m]<<"\t";
//               cout<<"t"<<temp<<"\t";
            }
//            cout<<"\n";
//            cout<<"t"<<temp<<"\t";
            index[i*(snpnum+1)+j]=index[i*(snpnum+1)+j]*anno_[j*anno_size+i];
            index[i*(snpnum+1)+snpnum]=index[i*(snpnum+1)+snpnum]-anno_[j*anno_size+i];
        }

//        cout<<"\n";
    }
}

struct param{
    int anno_size;
    int snpnum;
    int* anno_;
    double* index;
};

int Mstep_f (const gsl_vector* x, void* params, gsl_vector* f){

//    double *index=new double[anno_size*(snpnum+1)];
//    double *X=new double[anno_size];
//    Mstep_prepare(anno_size, post_prob, snpnum, k, anno_, index);


    int anno_size=((struct param *)params)->anno_size;
    int snpnum=((struct param *)params)->snpnum;
    int* anno_=((struct param *)params)->anno_;
    double* index=((struct param *)params)->index;

    double* X=new double[anno_size];
    double* Y=new double[snpnum];
    int i;
    int j;

    for(i=0;i<anno_size;i++){
        X[i]=gsl_vector_get(x,i);
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<anno_size;j++){
            Y[i]=0;
        }
    }

    for(i=0;i<snpnum;i++){
        for(j=0;j<anno_size;j++){
            Y[i]+=anno_[i*anno_size+j]*X[j];
        }
//        Y[i]=1/Y[i];
    }

    double* Y_=new double[snpnum+1];
    for(i=0;i<snpnum+1;i++){
        Y_[i]=1;
    }

    for(i=1;i<snpnum;i++){
        Y_[0]=Y_[0]*Y[i];
    }

    for(i=1;i<snpnum;i++){
        for(j=0;j<i;j++){
            Y_[i]=Y_[i]*Y[j];
        }
        for(j=i+1;j<snpnum;j++){
            Y_[i]=Y_[i]*Y[j];
        }
    }

    for(i=0;i<snpnum;i++){
        Y_[snpnum]=Y_[snpnum]*Y[i];
    }

//    delete [] Y;

    double* Z=new double[anno_size];

    for(i=0;i<anno_size;i++){
        for(j=0;j<=snpnum;j++){
            Z[i]=0;
        }
    }

    for(i=0;i<anno_size;i++){
        for(j=0;j<=snpnum;j++){
            Z[i]+=index[i*(snpnum+1)+j]*Y_[j];
        }
    }

//    delete [] Y_;

    for(i=0;i<anno_size;i++){
        gsl_vector_set(f,i,Z[i]);
    }

    return GSL_SUCCESS;
}

void Mstep_solve(int snpnum, int k, int anno_size, int* anno_, double* post_prob, double* new_gamma){

    double* index= new double[anno_size*(snpnum+1)];
    Mstep_prepare(anno_size, post_prob, snpnum, k, anno_, index);
//    int m;
//    int j;
/*    for(m=0;m<anno_size;m++){
        for(j=0;j<=snpnum;j++){
            cout<<index[m*(snpnum+1)+j]<<"\t";
        }
        cout<<"\n";
    }
*/
    const size_t n=anno_size;
    struct param p={anno_size,snpnum,anno_,index};
    gsl_multiroot_function f={&Mstep_f,n,&p};

    double *x_init=new double[anno_size];
    int i;
    for(i=0;i<anno_size;i++){
        x_init[i]=new_gamma[i];
    }

    gsl_vector *x=gsl_vector_alloc(n);
    for(i=0;i<anno_size;i++){
        gsl_vector_set(x,i,x_init[i]);
    }

    const gsl_multiroot_fsolver_type *T;
    T=gsl_multiroot_fsolver_hybrids;

    gsl_multiroot_fsolver *s;
    s=gsl_multiroot_fsolver_alloc(T,anno_size);

    gsl_multiroot_fsolver_set(s,&f,x);
    int status;
    int iter=0;
    //print_state(iter,s);
    cout<<iter<<"\t"<<s<<"\n";

    do{
        iter++;
        status=gsl_multiroot_fsolver_iterate(s);

        cout<<iter<<"\t";
        for(i=0;i<anno_size;i++){
            cout<<gsl_vector_get(s->x,i)<<"\t";
        }
        cout<<"\n";

        if(status)
            break;
        status=gsl_multiroot_test_residual(s->f,1e-7);
    }
    while(status==GSL_CONTINUE&&iter<1000);

    for(i=0;i<anno_size;i++){
        new_gamma[i]= gsl_vector_get(s->x,i);
    }

    printf("status=%s\n",gsl_strerror(status));

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
//    return 0;
}

