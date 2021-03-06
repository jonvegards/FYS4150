// Project 1 main program
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include "armadillo"    // How do I compile with this included?
#include "arma_solve.h"
#include "num_solve.h"
#include "lu_decomposition.h"

using namespace std;
using namespace arma;

void exact( double *, double *, int );
void error( double *, double *, double *, int );
void save_results( double *, double *, double *, double*, int );
void save_results_arma(mat, mat, int);
void error_lu( double *, double *, mat , int);

int main()
{
    int n = 100; // Program runs w/ n=10,100,1000,10000.
    double x[n+1], h, a[n+1], b[n+1], c[n+1], v[n+2], f[n+1], u_exact[n+1], err[n+1];
    int i;

    /////////////////////////////////////////////////////////////////////////////
    // Initializing arrays etc.
    /////////////////////////////////////////////////////////////////////////////

    h = 1 / (float(n)+1);
    // Imposing boundary conditions
    v[0] = 0; v[n] = 0; x[0] = 0; x[n+1] = 1;

    u_exact[n] = 0; err[n+1] = 0;
    a[0] = c[n+1] = 0;

    for (i=1; i <= n; i++){
        x[i] = i*h;
        f[i] = h*h*100*exp(-10*x[i]);
    }
    for (i=1; i <= n; i++){
        c[i] = -1.;
        a[i] = -1.;
        b[i] = 2.;
    }

    /////////////////////////////////////
    // Initializing matrices
    /////////////////////////////////////

    mat x_mat = zeros<mat>(n+1,1);
    mat x_mat2 = zeros<mat>(n+1,1);
    mat x2 = zeros<mat>(n+1,1);
    mat x3 = zeros<mat>(n+1,1);

    /////////////////////////////////////
    // Solving with our tridiagonal solver

    num_solve(n, a, b, c, v, f);

    /////////////////////////////////////
    // Solving with Armadillo

    for (i=1; i<=n; i++){
        x_mat(i) = i*h;
    }
    x_mat.reshape(n+2,1);
    x_mat(n+1) = 1.;

    x2 = arma_solve(n, h, x_mat, x2);

    /////////////////////////////////////
    // Solving by LU-decomposition

    x3 = lu_decomposition(n, h, x_mat, x3);

    /////////////////////////////////////

    // Finding the exact result
    exact( x, u_exact, n);

    // Computing the error
    error(err, u_exact, v, n);

    // Writing results to file
    //save_results( x, v, u_exact, err, n);
    //save_results_arma(x_mat, x2, n);
    return 0;
}

void exact( double *x, double *u_exact, int n){
    int i;
    // f(x) = 100*exp(-10*x)
    for (i=1; i<=n; i++){
        u_exact[i] = 1 - ( 1 - exp(-10) )*x[i] - exp(-10*x[i]);
    }
}

void error( double *err, double *u, double *v, int n){
    int i;
    double max;
    // Setting first element to 0.
    err[0] = 0;
    for (i=1; i<n; i++){
        // Calc. error using the given formula
        err[i] = log10( fabs( (v[i] - u[i]) / u[i] ) );
        if (fabs(err[i]) > fabs(err[i-1])){    // Picking out max error
            max = err[i];
        }
    }
    cout << "max error: " << max << endl;      // Printing max value
}

void save_results( double *x, double *v, double *u, double *err, int n){
    FILE *output_file;
    output_file = fopen("oppgave_b_n_100000.txt" , "w") ;  // Is there a way to produce several output files with different names?
    fprintf(output_file, "   %s    %s    %s    %s \n", "x", "v_numerical", "u", "error");
    int i;
    for (i=0; i<=n+1; i++){
        fprintf(output_file, "%12.5f %12.5f %12.5f %12.5f \n",
                 x[i], v[i], u[i], err[i] );
    }
    fclose (output_file);
}

void save_results_arma(mat err2, mat x2, int n){
    int i;
    ofstream myfile;
    myfile.open ("arm_solve_10000.txt");
    myfile << "x_mat" << "     " << "x2" << endl;
    for (i=0; i<n+1; i++){
        myfile << err2(i) << "    " << x2(i) << endl;
    }
    myfile.close();
}
