#include <iostream>
#include <time.h>
#include "armadillo"
#include "num_solve.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the problem with a tridiagonal
// solving algorithm
////////////////////////////////////////////////

mat num_solve(int n, mat a, mat b, mat c, mat v, mat f){
    // Copied from lecture notes p. 186
    // Note that the endpoints of f and v are sent into this
    // function, so we start indexing at 1 instead of 0.
    // First: forward substitution
    vec temp(n+1);

    double btemp = b(1);
    v(1) = f(1)/btemp;

    for(int i=1; i < n ; i++) {
        temp(i) = c(i-1)/btemp;
        btemp = b(i)-a(i)*temp(i);
        v(i) = (f(i) - a(i)*v(i-1))/btemp;
    }

    // Secondly: backsubstitution
    for(int i=n-1 ; i >= 1 ; i--) {
        v(i) -= temp(i+1)*v(i+1);
    }
    return v;
}
