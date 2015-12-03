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

void SkriveUtResultat(mat, int);

void SavingResultForTwoMoments(mat, mat, mat, int);
void SavingResultForTwoMoments2(mat, mat, mat, int);
void SavingResultForTwoMoments3(mat, mat, mat, int);

int main()
{
    int n = 99;
    double h = 1/((double) n + 1);
    double delta_t = h*h/2;
    double tsteps = 10000;
    double alpha = delta_t / (h*h);
    double alphaNC = (delta_t + delta_t/2) / (h*h); // Alpha-value for Crank-Nicolson
    int M = 2;                        // No. of samples of solution
    int sample1 = 100; int sample2 = 2990;
    double v_max = 0, v_0 = 1;        // Boundary conditions

    mat x = zeros<mat>(n+1,1);
    for (int i = 1; i <= n; i++) {
        x(i) = i*h;
    }

    //--------------------------------------------
    // Explicit scheme
    //--------------------------------------------
    mat v_explicit = zeros<mat>(n+1,M);
    mat v = zeros<mat>(n+1,1);
    v(0) = v_0;
    v(n) = v_max;

    // Time integration, for each time-iteration we obtain the new v-vector which we
    // can plot vs. position. We save the v-vector for each 100th round.
    int j = 0; // dummy variable for if-test
    for (int t = 0; t <= tsteps-1; t++) {
        for (int i = 1; i < n; i++) {
            // Discretized diff eq
            v(i) = alpha * v(i-1) + (1 - 2*alpha) * v(i) + alpha * v(i+1);
        }
        // If-test for saving results in matrix
        if(t == sample1 || t == sample2){
            v_explicit.col(j) = v;
            j++;
        }
    }


    SavingResultForTwoMoments(v_explicit.col(0), v_explicit.col(1), x, n);

    /*--------------------------------------------
    * Implicit scheme
    *
    * Initializing arrays etc.
    * Solving the eq mat(A)*vec(v) = vec(b) for vec(v)
    * where mat(A) has b on main diagonal, and a and c
    * on 1st upper/lower diagonal.
    *--------------------------------------------*/
    mat v_implicit = zeros<mat>(n+1,M);
    mat f = zeros<mat>(n+1,1);
    vec a(n+1), b(n+1), c(n+1), v_tri(n+1);

    // Imposing boundary conditions
    f(0) = v_0;
    f(n) = v_max;
    v_tri(0) = v_0;

    double b_val = 1+2*alpha;
    double val = -alpha;
    for (int i=1; i <= n; i++){
        c(i) = val;
        a(i) = val;
        b(i) = b_val;
    }
    a(0) = c(n) = a(n) = c(0) = 0;

    int q=0;
    for (int t = 0; t <= tsteps-1; t++) {
        v_tri = num_solve(n, a, b, c, v_tri, f);
        f(0) = v_0;
        f(n) = 0;
        f = v_tri;
        // If-test for saving results in matrix v_implicit
        if(t == sample1 || t == sample2){
            v_implicit.col(q) = v_tri;
            q++;
        }
    }

    SavingResultForTwoMoments2(v_implicit.col(0), v_implicit.col(1), x, n);

    /*--------------------------------------------
    * Crank-Nicolson scheme
    *
    * Initializing arrays etc.
    * Solving the eq mat(A_explicit)*vec(v_j) = mat(A_implicit)vec(v_{j-1}).
    * We solve the RHS by using the implicit schem, we
    * then obtain mat(A_explicit)*vec(v_j) = V_tilde,
    * which we solve by the explicit scheme.
    *--------------------------------------------*/
    mat v_tilde = zeros<mat>(n+1,1);
    mat v_CN = zeros<mat>(n+1,M);
    mat v_CNE = zeros<mat>(n+1,1);
    v_CNE(0) = v_0;
    v_tilde(0) = v_0;
    v_tilde(n) = v_max;

    b_val = 2+2*alphaNC;
    val = -alphaNC;
    for (int i=1; i <= n; i++){
        c(i) = val;
        a(i) = val;
        b(i) = b_val;
    }
    a(0) = c(n) = a(n) = c(0) = 0;

    // Time integration, for each time-iteration we obtain the new v-vector which we
    // can plot vs. position.
    q = 0; // dummy variable for if-test
    for (int t = 0; t <= tsteps-1; t++) {
        //v_tilde = ExplicitScheme(n, alpha, v_tilde);
        for (int i = 1; i < n; i++) {
            // Discretized diff eq
            v_CNE(i) = alphaNC * v_tilde(i-1) + (2 - 2*alphaNC) * v_tilde(i) + alphaNC * v_tilde(i+1);
        }
        v_tilde = v_CNE;
        v_tilde(0) = v_0;
        v_tilde(n) = v_max;
        v_CNE = num_solve(n, a, b, c, v_CNE, v_tilde);
        v_CNE(0) = v_0;
        v_CNE(n) = v_max;
        v_tilde = v_CNE;
        // If-test for saving results in matrix v_implicit
        if(t == sample1 || t == sample2){
            v_CN.col(q) = v_CNE;
            q++;
        }
    }

    SavingResultForTwoMoments3(v_CN.col(0), v_CN.col(1), x, n);

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
    output_file = fopen("oppgave_d_explicit.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
}

void SavingResultForTwoMoments2(mat V1, mat V2, mat x,int N){
    FILE *output_file;
    output_file = fopen("oppgave_d_implicit.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
}

void SavingResultForTwoMoments3(mat V1, mat V2, mat x,int N){
    FILE *output_file;
    output_file = fopen("oppgave_d_CN.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
}
