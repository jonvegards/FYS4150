// Program for solving the equation set by LU-decomposition
#include <iostream>
#include <cmath>


using namespace std;
//ludcmp(double a, int n, int indx, double &d);
//lubksb(double a, int n, int indx, double w);

int main(){
    /*int n = 1000, i, j, indx;
    double x[n], h, A[n][n], b[n], c[n], v[n], f[n], u_exact[n], err[n], d;

    h = 1 / (float(n)+1);
    v[0] = 0; v[n] = 0;

    // Creating matrix A

    for (i=1; i < n; i++){
        for (j=1; j < n; j++){
            if (fabs(j-i) < 1){
                A[i][j] = 2/(h*h);
            }
            if (fabs(j-i) == 1){
                A[i][j] = -1/(h*h);
            }
            else{
                A[i][j] = 0;
            }

        }
        x[i] = i*h;
        f[i] = 100*exp(-10*x[i]);
    }

    // LU-decomposition function decomposes A.
    //ludcmp( A, n, indx, &d);


    // Finding the solution x.
    //lubksb( A, n, indx, v);
*/
    return 0;
}
