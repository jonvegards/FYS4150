// Project 1 main program
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>    // For printing max value in an array
#include <armadillo>    // How do I compile with this included?
#include "lu_decomposition.h"

using namespace std;
using namespace arma;

void exact( double *, double *, int );
void error( double *, double *, double *, int );
void save_results( double *, double *, double *, double*, int );

int main()
{
    int n = 10, i, j;
    double x[n+1], h, a[n+1], b[n+1], c[n+1], v[n+2], f[n+1], u_exact[n+1], err[n+1];
    double btemp, temp[n+1];

    h = 1 / (float(n)+1);
    v[0] = 0; v[n] = 0; x[0] = 0; x[n+1] = 1;
    u_exact[n+1] = 0; err[n+1] = 0;

    for (i=1; i <= n; i++){
        x[i] = i*h;
        f[i] = h*h*100*exp(-10*x[i]);
    }
    for (i=1; i <= n; i++){
        c[i] = -1.;
        a[i] = -1.;
        b[i] = 2.;
    }
    a[0] = c[n+1] = 0;

    // Finding the exact result
    exact( x, u_exact, n);

    // f(x) = 100*exp(-10*x)
    // Copied from lecture notes p. 186
    // First: forward substitution

    btemp = b[1];
    v[1] = f[1]/btemp;

    for(i=2; i <= n ; i++) {
        temp[i] = c[i-1]/btemp;
        btemp = b[i]-a[i]*temp[i];
        v[i] = (f[i] - a[i]*v[i-1])/btemp;
    }

    // Secondly: backsubstitution
    for(i=n-1 ; i >= 1 ; i--) {
        v[i] -= temp[i+1]*v[i+1];
    }
    // Computing the error
    error(err, u_exact, v, n);

    // Writing results to file
    save_results( x, v, u_exact, err, n);
    cout << n << endl;

    // Solving the same problem by LU-decomposition

    mat A(n,n);
    vec f2(n);

    for (i=1; i < n; i++){
        f2(i) = h*h*100*exp(-10*x[i]);
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

    vec y = solve(A,f2);

    cout << y << endl;

    return 0;
}

void exact( double *x, double *u_exact, int n){
    int i;
    for (i=1; i<=n; i++){
        u_exact[i] = 1 - ( 1 - exp(-10) )*x[i] - exp(-10*x[i]);
    }
}

void error( double *err, double *u, double *v, int n){
    int i;
    double max;
    err[0] = 0;
    for (i=1; i<n; i++){
        err[i] = log10( fabs( (v[i] - u[i]) / u[i] ) );
        if (fabs(err[i]) > fabs(err[i-1])){
            max = err[i];
        }
    }
    cout << max << endl;                                    // Printing max value
}

void save_results( double *x, double *v, double *u, double *err, int n){
    FILE *output_file;
    output_file = fopen("oppgave_b_n_10.txt", "w") ;  // Is there a way to produce several output files with different names?
    fprintf(output_file, "   %s    %s    %s    %s \n", "x", "v", "u", "error");
    int i;
    for (i=0; i<=n+1; i++){
        fprintf(output_file, "%12.5f %12.5f %12.5f %12.5f \n",
                 x[i], v[i], u[i], err[i] );
    }
    fclose (output_file);
}
