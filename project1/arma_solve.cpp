#include <iostream>
#include <time.h>
#include "armadillo"
#include "arma_solve.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the same problem with Armadillo
////////////////////////////////////////////////

mat arma_solve(int n, double h, mat x_mat, mat x2){
    int i;

    clock_t start , finish ; // declare start and final time start = clock () ;

    // Declaring matrices
    mat A = zeros<mat>(n+1,n+1);
    mat f2 = zeros<mat>(n+1,1);

    // Initializing source term-vector
    for (i=1; i<=n; i++){
        f2(i) = h*h*100*exp(-10*x_mat(i));
    }

    // Filling matrix A
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1.;

    // Starting clock
    start = clock();

    // Calling on Armadillo's solve-function
    // to solve A * x2 = f2
    x2 = solve(A, f2);

    // Stopping clock
    finish = clock();

    printf ("Elapsed time Armadillo: %5.9f seconds.\n", ( (float( finish - start )) / CLOCKS_PER_SEC ));

    // Shifting the array one element in order to get the
    // boundary conditions correct
    x2.reshape(n+3,1);
    for (i=n; i>=0; i--){
        x2(i+1) = x2(i);
    }

    x2(0) = 0.;

    return x2;
}
