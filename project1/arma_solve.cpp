#include <iostream>
#include "armadillo"    // How do I compile with this included?
#include "arma_solve.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////////////
// Solving the same problem by LU-decomposition
/////////////////////////////////////////////////////////////////////////

double arma_solve(int n, double h, double *b, double *c){
    int i, j;

    mat A = zeros<mat>(n+1,n+1);
    mat f2 = zeros<mat>(n+1,1);
    mat x_mat = zeros<mat>(n+1,1);

    for (i=1; i<=n; i++){
        x_mat(i) = i*h;
    }

    for (i=1; i <= n; i++){
        f2(i) = h*h*100*exp(-10*x_mat(i));
        for (j=1; j<=n; j++){
            if (fabs(i-j) == 1){
                A(i,j) = c[i];
            }
            if ( (i-j) < 0 + 1.0e-1 and (i-j) > 0 - 1.0e-1){
                A(i,j) = b[i];
            }
        }
    }
    A(0,0) = 2.;
    A(0,1) = -1.;
    A(1,0) = -1.;

    //cout << A << endl;

    mat x2 = solve(A, f2);

    x2.reshape(n+3,1);
    x_mat.reshape(n+3,1);
    x_mat(n+1) = 1.;

    // Shifting the array one element
    for (i=n; i>0; i--){
        x2(i) = x2(i-1);
    }

    x2(0) = 0.;

    ofstream myfile;
    myfile.open ("arm_solve_1000.txt");
    myfile << "x_mat" << "     " << "x2" << endl;
    for (i=0; i<=n+1; i++){
        myfile << x_mat(i) << "    " << x2(i) << endl;
    }
    myfile.close();

    cout << "Success" << endl;

    return 0;
}
