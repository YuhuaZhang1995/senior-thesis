#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

double mvnnorm(int* n,
               double* lower,
               double* upper,
               int* infin,
               double* correl,
               int* maxpts,    // param
               double* abseps, // param
               double* releps, // param
               double* error,  // estimated abs. error. with 99% confidence interval
               double* value,     // results store here.
               int* inform)    // inform message goes here
{
    mvndst_(n,lower, upper, infin, correl,
          maxpts, abseps, releps, error, value, inform);
//    printf ("error = %g, value = %g, inform = %d\n", *error, *value, *inform);
    return *value;
}

double pmvnorm(int n,
                 double* bound,
                 double* correlationMatrix) // (2,1), (3,1), (3,2) .....
{
  int maxpts_ = 25000;     // default in mvtnorm: 25000
  double abseps_ = 1e-6;   // default in mvtnorm: 0.001, we make it more stringent
  double releps_ = 0;      // default in mvtnorm: 0

  double* lower = new double[n];
  int* infin = new int[n];

  int i = 0;
  for (i = 0; i < n; ++i) {
    infin[i] = 0; // (-inf, bound]
    lower[i] = 0;
  }

  // return values
  double value_ = 0;
  double error=0;
  int inform_ = 0;

  double ret = mvnnorm(&n, lower, bound, infin, correlationMatrix, &maxpts_, &abseps_, &releps_, &error, &value_, &inform_);
  delete[] (lower);
  delete[] (infin);

  return ret;
}


