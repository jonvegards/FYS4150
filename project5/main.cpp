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
#include "lib.h"

using namespace std;
using namespace arma;

void SavingResultForTwoMoments(mat, mat, mat, int);
void SavingResultForTwoMoments2(mat, mat, mat, int);
void SavingResultForTwoMoments3(mat, mat, mat, int);
// function for gaussian random numbers
double gaussian_deviate(long *);

void MC_simulation(double tsteps, int max_trials, double delta_t, double move_probability, int n);
void MC_simulation_gaussian(double tsteps, int max_trials, double delta_t, double move_probability, int n);

int main()
{
    int n = 99;
    double h = 1/((double) n + 1);
    double delta_t = h*h/2;
//    double tsteps = 5000;
//    double alpha = delta_t / (h*h);
//    double alphaNC = (delta_t + delta_t/2) / (h*h); // Alpha-value for Crank-Nicolson
//    int M = 2;                        // No. of samples of solution
//    int sample1 = 200; int sample2 = 5000;
//    double v_max = 0, v_0 = 1;        // Boundary conditions

//    mat x = zeros<mat>(n+1,1);
//    for (int i = 1; i <= n; i++) {
//        x(i) = i*h;
//    }

//    //--------------------------------------------
//    // Explicit scheme
//    //--------------------------------------------
//    mat v_explicit = zeros<mat>(n+1,M);
//    mat v = zeros<mat>(n+1,1);
//    v(0) = v_0;
//    v(n) = v_max;

//    // Time integration, for each time-iteration we obtain the new v-vector which we
//    // can plot vs. position. We save the v-vector for each 100th round.
//    int j = 0; // dummy variable for if-test
//    for (int t = 1; t <= tsteps; t++) {
//        for (int i = 1; i < n; i++) {
//            // Discretized diff eq
//            v(i) = alpha * v(i-1) + (1 - 2*alpha) * v(i) + alpha * v(i+1);
//        }
//        // If-test for saving results in matrix
//        if(t == sample1 || t == sample2){
//            v_explicit.col(j) = v;
//            j++;
//        }
//    }


//    SavingResultForTwoMoments(v_explicit.col(0), v_explicit.col(1), x, n);

//    /*--------------------------------------------
//        * Implicit scheme
//        *
//        * Initializing arrays etc.
//        * Solving the eq mat(A)*vec(v) = vec(b) for vec(v)
//        * where mat(A) has b on main diagonal, and a and c
//        * on 1st upper/lower diagonal.
//        *--------------------------------------------*/
//    mat v_implicit = zeros<mat>(n+1,M);
//    mat f = zeros<mat>(n+1,1);
//    vec a(n+1), b(n+1), c(n+1), v_tri(n+1);

//    // Imposing boundary conditions
//    f(0) = v_0;
//    f(n) = v_max;
//    v_tri(0) = v_0;

//    double b_val = 1+2*alpha;
//    double val = -alpha;
//    for (int i=1; i <= n; i++){
//        c(i) = val;
//        a(i) = val;
//        b(i) = b_val;
//    }
//    a(0) = c(n) = a(n) = c(0) = 0;

//    int q=0;
//    for (int t = 1; t <= tsteps; t++) {
//        v_tri = num_solve(n, a, b, c, v_tri, f);
//        f(0) = v_0;
//        f(n) = 0;
//        f = v_tri;
//        // If-test for saving results in matrix v_implicit
//        if(t == sample1 || t == sample2){
//            v_implicit.col(q) = v_tri;
//            q++;
//        }
//    }

//    SavingResultForTwoMoments2(v_implicit.col(0), v_implicit.col(1), x, n);

//    /*--------------------------------------------
//        * Crank-Nicolson scheme
//        *
//        * Initializing arrays etc.
//        * Solving the eq mat(A_explicit)*vec(v_j) = mat(A_implicit)vec(v_{j-1}).
//        * We solve the RHS by using the implicit schem, we
//        * then obtain mat(A_explicit)*vec(v_j) = V_tilde,
//        * which we solve by the explicit scheme.
//        *--------------------------------------------*/
//    mat v_tilde = zeros<mat>(n+1,1);
//    mat v_CN = zeros<mat>(n+1,M);
//    mat v_CNE = zeros<mat>(n+1,1);
//    v_CNE(0) = v_0;
//    v_tilde(0) = v_0;
//    v_tilde(n) = v_max;

//    b_val = 2+2*alphaNC;
//    val = -alphaNC;
//    for (int i=1; i <= n; i++){
//        c(i) = val;
//        a(i) = val;
//        b(i) = b_val;
//    }
//    a(0) = c(n) = a(n) = c(0) = 0;

//    // Time integration, for each time-iteration we obtain the new v-vector which we
//    // can plot vs. position.
//    q = 0; // dummy variable for if-test
//    for (int t = 1; t <= tsteps; t++) {
//        //v_tilde = ExplicitScheme(n, alpha, v_tilde);
//        for (int i = 1; i < n; i++) {
//            // Discretized diff eq
//            v_CNE(i) = alphaNC * v_tilde(i-1) + (2 - 2*alphaNC) * v_tilde(i) + alphaNC * v_tilde(i+1);
//        }
//        v_tilde = v_CNE;
//        v_tilde(0) = v_0;
//        v_tilde(n) = v_max;
//        v_CNE = num_solve(n, a, b, c, v_CNE, v_tilde);
//        v_CNE(0) = v_0;
//        v_CNE(n) = v_max;
//        v_tilde = v_CNE;
//        // If-test for saving results in matrix v_implicit
//        if(t == sample1 || t == sample2){
//            v_CN.col(q) = v_CNE;
//            q++;
//        }
//    }

//    SavingResultForTwoMoments3(v_CN.col(0), v_CN.col(1), x, n);

    /* Algorithm for Monte Carlo-simulation
     *
     * // Declare a matrix whose column vectors describe the system,
     * // i.e. the first column has N particles always, then the
     * // next column contain particles that has diffused one step
     * // and so on.
     * }
     */

    int max_trials= 1e6;
    int tsteps = 200;
    double move_probability = 0.5;
    //MC_simulation(tsteps, max_trials, delta_t, move_probability, n);
    MC_simulation_gaussian(tsteps, max_trials, delta_t, move_probability, n);

    return 0;
} // End main function

void SavingResultForTwoMoments(mat V1, mat V2, mat x,int N){
    // Function for writing results from explicit method
    FILE *output_file;
    output_file = fopen("oppgave_d_explicit.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
} // End SavingResultForTwoMoments

void SavingResultForTwoMoments2(mat V1, mat V2, mat x,int N){
    // Function for writing results from implicit method
    FILE *output_file;
    output_file = fopen("oppgave_d_implicit.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
} // End SavingResultForTwoMoments2

void SavingResultForTwoMoments3(mat V1, mat V2, mat x,int N){
    // Function for writing results from Crank-Nicolson method
    FILE *output_file;
    output_file = fopen("oppgave_d_CN.txt" , "w") ;  // Is there a way to produce several output files with different names?

    fprintf(output_file, "   %s    %s    %s\n", "x", "V1", "V2");
    for (int i=0; i<=N; i++){
        fprintf(output_file, "%f %f %f \n",
                x(i), V1(i), V2(i) );
    }
    fclose (output_file);
} // End SavingResultForTwoMoments3

// random numbers with gaussian distribution
double gaussian_deviate(long * idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran2(idum) -1.0;
      v2 = 2.*ran2(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
} // end function for gaussian deviates

void MC_simulation(double tsteps, int max_trials, double delta_t, double move_probability, int n){
    double steplength = sqrt(2*delta_t);
    mat walk_cumulative = zeros<mat>(max_trials+1,1);
    mat walk2_cumulative = zeros<mat>(max_trials+1,1);
    mat probability_vec = zeros<mat>(n+1,1);
    probability_vec[0] = max_trials;

    long idum = -1;                         // initialising random number generator
    for (int trial=1; trial <= max_trials; trial++){
        int position = 0;
        for (int walks = 1; walks <= tsteps; walks++){
            if (ran1(&idum) <= move_probability){
                position += 1;
            }
            else{
                position -= 1;
            }
            if(position <= 0 or position > 99) break;
            walk2_cumulative[trial] += position;//*position;
            probability_vec[position] += 1;
        } // end of loop over walks
    } // end of loop over trials

    cout << "total time:" << delta_t*tsteps << endl;

    FILE *file_walk;
    file_walk = fopen("oppgave_MC_walk200.txt" , "w") ;
    fprintf(file_walk, "%s       %s       %s      \n", "i","position","variance");
    for(int i=1; i<=max_trials; i++){
        double xaverage = walk_cumulative[i]/((double) max_trials);
        double xaverage2 = walk2_cumulative[i]/((double) max_trials*max_trials);
        double variance = xaverage2 - xaverage*xaverage;
        fprintf(file_walk, "%d %.8f %.8f %.8f \n", i,xaverage,variance, walk_cumulative[i]/(double) max_trials);
        // Print to file
    }
    fclose (file_walk);

    FILE *output_file;
    output_file = fopen("oppgave_MC_probability200.txt" , "w") ;
    fprintf(output_file, "   %s  \n", "histogram");
    for(int i=0; i<n; i++){
        double histogram = probability_vec[i]/probability_vec[0];
        fprintf(output_file, "%d %.8e\n",i, histogram);
        // Print to file
    }
    fclose (output_file);
} // end function MC_simulation

void MC_simulation_gaussian(double tsteps, int max_trials, double delta_t, double move_probability, int n){
    long idum2 = -1;
    double steplength=0;
    mat walk_cumulative = zeros<mat>(max_trials+1,1);
    mat walk2_cumulative = zeros<mat>(max_trials+1,1);
    mat probability_vec = zeros<mat>(n+1,1);
    probability_vec[0] = max_trials;
    double h = 1/((double) n + 1); // width of x-slice

    /*
     * We have divided the interval [0,1] into 100 slices,
     * so we have to check if the steplength is big enough
     * to move the particle into the neighbouring slice
     */

    long idum = -1;                         // initialising random number generator
    for (int trial=1; trial <= max_trials; trial++){
        int position = 0;
        for (int walks = 1; walks <= tsteps; walks++){
            steplength = sqrt(2*delta_t)*gaussian_deviate(&idum2);
            if(ran1(&idum) <= move_probability){
                // Move right
                if(steplength > h) position += 1;
            }
            else{
                if(steplength > h) position -= 1;
            }
            if(position <= 0 or position > 99) break;
            walk2_cumulative[trial] += position;//*position;
            probability_vec[position] += 1;
        } // end of loop over walks
    } // end of loop over trials

    cout << "total time:" << delta_t*tsteps << endl;

    FILE *file_walk;
    file_walk = fopen("oppgave_MC_walk200gauss.txt" , "w") ;
    fprintf(file_walk, "%s       %s       %s      \n", "i","position","variance");
    for(int i=1; i<=max_trials; i++){
        double xaverage = walk_cumulative[i]/((double) max_trials);
        double xaverage2 = walk2_cumulative[i]/((double) max_trials*max_trials);
        double variance = xaverage2 - xaverage*xaverage;
        fprintf(file_walk, "%d %.8f %.8f %.8f \n", i,xaverage,variance, walk_cumulative[i]/(double) max_trials);
        // Print to file
    }
    fclose (file_walk);

    FILE *output_file;
    output_file = fopen("oppgave_MC_probability200gauss.txt" , "w") ;
    fprintf(output_file, "   %s  \n", "histogram");
    for(int i=0; i<n; i++){
        double histogram = probability_vec[i]/probability_vec[0];
        fprintf(output_file, "%d %.8e\n",i, histogram);
        // Print to file
    }
    fclose (output_file);
} // End MC_simulation_gaussian function
