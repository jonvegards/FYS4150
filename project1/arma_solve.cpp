#include <iostream>
#include "armadillo"    // How do I compile with this included?
#include "arma_solve.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////////////
// Solving the same problem by LU-decomposition
/////////////////////////////////////////////////////////////////////////

mat arma_solve(int n, double h, mat x_mat, mat x2){
    int i;

    mat A = zeros<mat>(n+1,n+1);
    mat f2 = zeros<mat>(n+1,1);

    for (i=1; i<=n; i++){
        f2(i) = h*h*100*exp(-10*x_mat(i));
    }

    // Filling matrix A
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1.;

    //cout << A << endl;

    x2 = solve(A, f2);

    x2.reshape(n+3,1);

    // Shifting the array one element
    for (i=n; i>=0; i--){
        x2(i+1) = x2(i);
    }

    x2(0) = 0.;

    cout << "Success" << endl;
    return x2;
}
