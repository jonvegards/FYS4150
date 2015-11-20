#include <iostream>
#include <time.h>
#include "armadillo"
#include "arma_solve.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the same problem with Armadillo
////////////////////////////////////////////////

mat arma_solve(int n, double h, mat x2, vec f2, double alpha){
    int i;

    clock_t start , finish ; // declare start and final time start = clock () ;

    // Declaring matrices
    mat A = zeros<mat>(n+1,n+1);
    mat B = zeros<mat>(n+1,n+1);
    mat identity = eye<mat>(n+1,n+1);

    // Filling matrix B
    B.diag() += 2.;
    B.diag(1) += -1.;
    B.diag(-1) += -1.;

    // Setting matrix A = identity + alpha * B
    A = identity + alpha*B;
    cout << A << endl;

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

    x2(0) = 1.;

    return x2;
}
