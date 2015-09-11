#include <iostream>
#include <time.h>
#include "num_solve.h"

using namespace std;

////////////////////////////////////////////////
// Solving the problem with a tridiagonal
// solving algorithm
////////////////////////////////////////////////

void num_solve(int n, double *a, double *b, double *c, double *v, double *f){
    int i;

    // Copied from lecture notes p. 186
    // First: forward substitution
    clock_t start , finish ; // declare start and final time start = clock () ;
    double temp[n+1];

    start = clock(); // Starting the timer

    double btemp = b[1];
    v[1] = f[1]/btemp;

    for(i=2; i <= n ; i++) {
        temp[i] = c[i-1]/btemp;
        btemp = b[i]-c[i]*temp[i];
        v[i] = (f[i] - a[i]*v[i-1])/btemp;
    }

    // Secondly: backsubstitution
    for(i=n-1 ; i >= 1 ; i--) {
        v[i] -= temp[i+1]*v[i+1];
    }
    finish = clock(); // Stopping the timer

    printf ("Elapsed time numerical: %5.9f seconds.\n", ( (float)( finish - start ) / CLOCKS_PER_SEC ));
}
