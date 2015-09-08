// Function for finding the LU decomposed matrix A

#include "lu_decomposition.h"

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <armadillo>

using namespace std;
using namespace arma;
/*
double lu_decomposition(double *a, double *b, double *c, double *x, int n){
    int i, j;
    mat A(n,n);
    vec f(n);
    double h = 1/(n+1);

    for (i=1; i <= n; i++){
        c[i] = -1.;
        a[i] = -1.;
        b[i] = 2.;
        f(i) = h*h*100*exp(-10*x[i]);
        for (j=1; j<n; j++){
                if (i==j){
                    A(i,j) = b[i];
                }
                if (i-j == 1){
                    A(i,j) = c[i];
                }
                else{
                    A(i,j) = 0.;
                }
        }
    }

        // Playing zone for Armadillo LU
        //mat P,L,U;
        //lu(L, U, P, A);

        vec y = solve(A,f);

        // And how do I use the arma-package?
    return 0;
}
*/
