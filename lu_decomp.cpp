#include <iostream>
#include <time.h>
#include "armadillo"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the same problem by LU-decomposition
////////////////////////////////////////////////

int main(int argc, char const *argv[]){
    int n=atof(argv[1]), i;
    mat w, y, L_inv, U_inv, x;
    double h = 1 / (float(n)+1);;

    clock_t start , finish ; // declare start and final time start = clock () ;
    start = clock();

    mat A = zeros<mat>(n+1,n+1);

    // Filling matrix A
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1.;

    mat f2 = zeros<mat>(n+1,1);
    mat x_mat = zeros<mat>(n+1,1);

    for (i=1; i<=n; i++){
        x_mat(i) = i*h;
    }
    x_mat.reshape(n+2,1);
    x_mat(n+1) = 1.;

    for (i=1; i<=n; i++){
        f2(i) = h*h*100*exp(-10*x_mat(i));
    }

    w = f2;

    mat L,U,P;
    lu(L,U,P,A);

    L_inv = inv(L);

    y = L_inv*w;

    U_inv = inv(U);

    x = U_inv * y;

    finish = clock();

    //x.print();

    printf ("Elapsed time numerical: %5.9f seconds.\n", ( (float)( finish - start ) / CLOCKS_PER_SEC ));

    cout << "Success" << endl;
    return 0;
}
