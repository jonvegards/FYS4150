#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the same problem by LU-decomposition
////////////////////////////////////////////////

int lu_decomp(){
    int n=1000;

    mat A = zeros<mat>(n+1,n+1);

    // Filling matrix A
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1.;

    mat L,U,P;
    lu(L,U,P,A);

    cout << "Success" << endl;
    return 0;
}
