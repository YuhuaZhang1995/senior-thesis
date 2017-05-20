#include "header.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

using namespace std;

int main(){
//    int samplesize=10;
//    int snpnum=5;
//    double X[]={0,1,2,2,1,1,2,1,2,0,0,0,0,1,2,0,2,1,2,1,2,0,0,0,2,0,0,0,1,2,2,1,0,1,0,0,0,1,2,0,1,0,2,2,0,0,1,0,0,0};
    int samplesize=2054;
    FILE *fin;
	fin=fopen("debug1","r");

	int snpnum=0;
	int c;
	for( ;; ){
		c=fgetc(fin);
		if(c=='\n')
			break;
		if(c=='.'){
			++snpnum;
		}
	}
//	snpnum=snpnum+1;
	cout<<snpnum<<"\n";

        fclose(fin);

	double *LD_mat=new double[snpnum*snpnum];
	ifstream f("debug1",std::ios_base::in);
	int i;
	int j;
	if(f.is_open()){
		cout<<"File opened success"<<"\n";
		for(i=0;i<snpnum*snpnum;i++){
			f>>LD_mat[i];
		}
	}
//    FILE* fin;
	fin=fopen("debug2","r");
	int count_=0;
	for( ;; ){
		c=fgetc(fin);
		if(c=='\t'){
			++count_;
		}
		if(c=='\n')
            break;
	}
	count_=count_+1;
//	cout<<count_<<'\n';
	fclose(fin);

	double *new_LD=new double[count_*count_];
	int *cnt=new int[count_];
	ifstream ff("debug2",std::ios_base::in);
	for(i=0;i<count_;i++){
		ff>>cnt[i];
	}
	int cnti;
	int cntj;

	for(i=0;i<count_;i++){
		for(j=0;j<count_;j++){
            cnti=cnt[i];
            cntj=cnt[j];
			new_LD[i*count_+j]=LD_mat[(cnti-1)*snpnum+cntj-1];
//			cout<<new_LD[i*count_+j]<<"\t";
		}
//		cout<<"\n";
	}

	snpnum=count_;
	count_=0;

	ifstream fff("zscore",std::ios_base::in);
	double *z_score=new double[snpnum];
	for(i=0;i<snpnum;i++){
		fff>>z_score[i];
//		cout<<z_score[i]<<"\t";
	}
//	cout<<"\n";

	fin=fopen("anno","r");
	for( ;; ){
		c=fgetc(fin);
		if(c=='\n')
			break;
		if(c!='\t'){
            ++count_;
		}
	}
//	cout<<count_<<"\n";
	fclose(fin);
	int s_anno=count_;
	count_=0;

	ifstream ffff("anno",std::ios_base::in);
	int *anno=new int[snpnum*s_anno];
	for(i=0;i<snpnum;i++){
		for(j=0;j<s_anno;j++){
			ffff>>anno[i*s_anno+j];
//			cout<<anno[i*s_anno+j]<<"\t";
		}
//		cout<<"\n";
	}

	int k=1;
	double* post_prob=new double[snpnum*(k+1)];
	double* old_prob=new double[snpnum];
	for(i=0;i<snpnum;i++){
		old_prob[i]=0;
	}
	double* mar_prob=new double[snpnum];
	double* gamma=new double[s_anno];
	for(i=0;i<s_anno;i++){
		gamma[i]=0.5;
	}
	int iter_max=10;

	for(i=0;i<iter_max;i++){
//		Estep_search(s_anno,mar_prob,samplesize,k,snpnum,z_score,new_LD,gamma,anno,post_prob);
        MCMC_smaple(mar_prob, post_prob, z_score, samplesize, new_LD, k, snpnum, gamma, anno, s_anno);
//        for(j=0;j<snpnum;j++){
//            cout<<mar_prob[j]<<"\t";
//        }
//        cout<<"\n";
		Mstep_solve(snpnum,k,s_anno,anno,post_prob,gamma);
		for(j=0;j<snpnum;j++){
			if(mar_prob[j]-old_prob[j]<0.01){
				count_++;
			}
		}
		if(count_==snpnum)
			break;
		count_=0;
		for(j=0;j<snpnum;j++){
			old_prob[j]=mar_prob[j];
		}
	}

	ofstream mydata;
	mydata.open("result.txt");
    for(i=0;i<snpnum;i++){
        mydata<<mar_prob[i]<<"\n";
    }
    mydata.close();






//	for(i=0;i<10;i++){
//		cout<<LD_mat[i]<<"\n";
//	}

//	double *sigma_k=new double[count];
//	for(i=0;i<count;i++){
//        sigma_k[i]=1;
//	}

//    double* CorrelationMatrix=new double[count*(count-1)/2];
//    double sigma_e=1;
//    cal_corr(new_LD, count, sigma_k, CorrelationMatrix,sigma_e,samplesize);

//    for(i=1;i<4;i++){
//        for(j=0;j<i;j++){
//            cout<<CorrelationMatrix[i*(i-1)/2+j]<<"\t";
//        }
//        cout<<"\n";
//    }

	return 0;


//    double* sigma=new double[snpnum*snpnum];
//    LD_mat(sigma, X, snpnum, samplesize);



//    double sigma_k[]={1,2,1,0,1};
//    double* CorrelationMatrix=new double[snpnum*(snpnum-1)/2];
//    double sigma_e=1;
//    cal_corr(sigma, snpnum, sigma_k, CorrelationMatrix,sigma_e,samplesize);

 /*   int i;
    int j;
    for(i=0;i<snpnum;i++){
        for(j=0;j<snpnum;j++){
            cout<<CorrelationMatrix[i*snpnum+j]<<"\t";
        }
        cout<<"\n";
    }*/
}
