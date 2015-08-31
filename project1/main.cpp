// Project 1 main program
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>    // For printing max value in an array
#include <armadillo>

using namespace std;

void exact( double *, double *, int );
void error( double *, double *, double *, int );
void save_results( double *, double *, double *, double*, int );

int main()
{
    int n = 1000, i, j;
    double x[n], h, a[n], b[n], c[n], v[n], f[n], u_exact[n], err[n];
    arma::mat A;
    double btemp, temp[n];

    h = 1 / (float(n)+1);
    v[0] = 0; v[n] = 0;

    for (i=1; i < n; i++){
        x[i] = i*h;
        f[i] = h*h*100*exp(-10*x[i]);
        c[i] = -1.;
        a[i] = -1.;
        b[i] = 2.;
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
    a[0] = c[n] = 0;

    // Finding the exact result
    exact( x, u_exact, n);

    // f(x) = 100*exp(-10*x)
    // Copied from lecture notes p. 186
    // First: forward substitution

    btemp = b[1];
    v[0] = f[0]/btemp;

    for(i=1; i < n ; i++) {
        temp[i] = c[i-1]/btemp;
        btemp = b[i]-a[i]*temp[i];
        v[i] = (f[i] - a[i]*v[i-1])/btemp;
    }

    // Secondly: backsubstitution
    for(i=n-1 ; i > 1 ; i--) {
        v[i] -= temp[i+1]*v[i+1];
    }
    // Computing the error
    error(err, u_exact, v, n);

    // Writing results to file
    save_results( x, v, u_exact, err, n);

    // Playing zone for Armadillo LU
    arma::mat P,L,U;
    arma::lu(L, U, P, A);

    return 0;
}

void exact( double *x, double *u_exact, int n){
    int i;
    for (i=1; i<n; i++){
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
    //double max = *max_element(err, err + n);            // This doesn't work properly.
    //cout << max << endl;
}

void save_results( double *x, double *v, double *u, double *err, int n){
    FILE *output_file;
    output_file = fopen("oppgave_b_n_1000.txt", "w") ;  // Is there a way to produce several output files with different names?
    fprintf(output_file, "%s %s %s %s \n", "x", "v", "u", "error");
    int i;
    for (i=0; i<n; i++){
        fprintf(output_file, "%12.5f %12.5f %12.5f %12.5f \n",
                 x[i], v[i], u[i], err[i] );
    }
    fclose (output_file);
}
