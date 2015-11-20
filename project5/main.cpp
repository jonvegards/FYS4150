/* Project 5, FYS4150
 * Reusing code from project 1
 *
 * Solving the 1D-diffusion equation
 * with explicit and implisit Euler-methods,
 * Crank-Nicolson and finally simulating
 * with Monte Carlo-methods.
 */

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include "armadillo"    // How do I compile with this included?
#include "arma_solve.h"
#include "num_solve.h"

using namespace std;
using namespace arma;

void exact( double *, double *, int );
void error( double *, double *, double *, int );
void save_results( double *, double *, double *, double*, int );
void save_results_arma(mat, mat, int);
void error_lu( double *, double *, mat , int);

void SkriveUtResultat(mat, int);
void SavingResultForTwoMoments(mat, mat, mat, int);

int main()
{
    int n = 10; // Program runs w/ n=10,100,1000,10000.

    //--------------------------------------------
    // Implicit scheme
    //--------------------------------------------
    double delta_t = (h*h)/4;
    double tsteps = 50;
    double alpha = delta_t / (h*h);
    int M = 2;                        // No. of samples of solution
    mat v = zeros<mat>(n+2,1);
    mat vnew = zeros<mat>(n+2,M);
    mat x = zeros<mat>(n+2,1);

    v(n+1) = vnew(n+1) = 0.0;
    v(0) = 1;
    x(n+1) = (n+1)*h;
    for (int i = 1; i < n+1; i++) {
        x(i) = i*h;
        // initial condition
        v(i) = 0;
        vnew(i) = 0;
        // intitialise the new vector unew(i) = 0;
    }

    // Time integration, for each time-iteration we obtain the new v-vector which we
    // can plot vs. position. We save the v-vector for each 100th round.
    int j = 0; // dummy variable for if-test
    for (int t = 0; t <= tsteps-1; t++) {
        for (int i = 1; i < n+1; i++) {
            // Discretized diff eq
            v(i) = alpha * v(i-1) + (1 - 2*alpha) * v(i) + alpha * v(i+1);
        }
        // If-test for saving results in matrix
        if(t == 1 || t == tsteps-1){
            for(int k=0;k<=n+1;k++){
                vnew(k,j) = v(k);
            }
            j++;
        }

    }

    SavingResultForTwoMoments(vnew.col(0), vnew.col(1), x, n+1);

    //--------------------------------------------
    // Explicit scheme
    //--------------------------------------------
    double x[n+1], h, a[n+1], b[n+1], c[n+1], v[n+2], f[n+1], u_exact[n+1], err[n+1];

    /////////////////////////////////////////////////////////////////////////////
    // Initializing arrays etc.
    // Solving the eq mat(A)*vec(v) = vec(b), where mat(A) has b on main
    // diagonal, and a and c on 1st upper/lower diagonal. f_i = h*h*f(x_i)
    /////////////////////////////////////////////////////////////////////////////

    // Imposing boundary conditions
    v[0] = 1; v[n] = 0; x[0] = 0; x[n+1] = 1;

    u_exact[n] = 0; err[n+1] = 0;
    a[0] = c[n+1] = 0;

    for (i=1; i <= n; i++){
        x[i] = i*h;
        f[i] = 0;
    }
    for (i=1; i <= n; i++){
        c[i] = -1.;
        a[i] = -1.;
        b[i] = 2.;
    }

    // Solving with our tridiagonal solver
    num_solve(n, a, b, c, v, f);

    /////////////////////////////////////
    // Initializing matrices
    /////////////////////////////////////

    //    mat x_mat = zeros<mat>(n+1,1);
    //    mat x_mat2 = zeros<mat>(n+1,1);
    //    mat x2 = zeros<mat>(n+1,1);
    //    mat x3 = zeros<mat>(n+1,1);

    /////////////////////////////////////
    // Solving with Armadillo

    //    for (i=1; i<=n; i++){
    //        x_mat(i) = i*h;
    //    }
    //    x_mat.reshape(n+2,1);
    //    x_mat(n+1) = 1.;

    //    // Initializing source term-vector
    //    mat V = zeros<mat>(n+1,1);
    //    for (i=1; i<=n; i++){
    //        V(i) = 0;
    //    }

    //    x2 = arma_solve(n, h, x2, V, alpha);

    /////////////////////////////////////

    // Finding the exact result
    //exact( x, u_exact, n);

    // Computing the error
    //error(err, u_exact, v, n);

    // Writing results to file
    //save_results( x, v, u_exact, err, n);
    //save_results_arma(v, x, n);
    return 0;
}

void SkriveUtResultat(mat probability, int N){
    cout << "probability " << "= array([";
    for (int i=0; i<N; i++){
        cout << probability(i) << ", ";
    }
    cout << "])" << endl;
}

void SavingResultForTwoMoments(mat V1, mat V2, mat x,int N){
    FILE *output_file;
    output_file = fopen("oppgave_d.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
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
