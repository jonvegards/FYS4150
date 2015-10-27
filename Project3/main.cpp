// Program for calculating a 6D-integral
// with the methods Gauss-Legendre, Gauss-Laguerre,
// Monte Carlo (both brute force and not-so-
// brute-force). In the end we may parallelize
// the code

#include <iostream>
#include <cmath>
#include "lib.h"
#include "gauss-laguerre.cpp"
#include "time.h"

using namespace std;


void StupidIntegrationMethod_aka_Gauss_Legendre(int);
void ABitBetterIntegrationMethod_aka_Gauss_Laguerre(int);
void MonteCarloBruteForce(int n);
void MonteCarlo(int n);

void gauss_laguerre(double *x, double *w, int n, double alf);

double IntegrandSphericalCoordinates(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
double IntegrandMC(double *);
double Integrand(double x1, double y1, double z1, double x2, double y2, double z2);
double IntegrandMC_radial(double *);


// This program use exampleprogram.cpp as a template
int main(){
    // main() calls on other functions in order to use
    // the different integration methods; Gauss-Legendre,
    // Gauss-Laguerre, brute force Monte Carlo and a
    // improved Monte Carlo-method.

    double pi = 4*atan(1);
    clock_t start, finish;

    int *n_gauss;
    // Def. a list which contains no. of mesh points
    // It's wise to choose sth. below 30, if not the program
    // will be extremely slow.
    int M=5;
    n_gauss = new int[M];
    n_gauss[0] = 10; n_gauss[1]=15; n_gauss[2]=20; n_gauss[3] = 25; n_gauss[4] = 30;
    printf("Gauss-Legendre:\n");
    printf("Integral:   Relative error:     n: \n");
    for(int i=0; i<M; i++){
        start = clock();
        StupidIntegrationMethod_aka_Gauss_Legendre(n_gauss[i]);
        finish = clock();
        cout << "Elapsed time Gauss-Legendre: " << ((double) (finish-start)/CLOCKS_PER_SEC) << "  for n=" << n_gauss[i] << endl;
    }
    printf("Gauss-Laguerre:\n");
    printf("Integral:   Relative error:     n: \n");
    for(int i=0; i<M; i++){
        start = clock();
        ABitBetterIntegrationMethod_aka_Gauss_Laguerre(n_gauss[i]);
        finish = clock();
        cout << "Elapsed time Gauss-Laguerre: " << ((double) (finish-start)/CLOCKS_PER_SEC) << "  for n=" << n_gauss[i] << endl;
    }

    // Part for Monte Carlo-methods
    int *n_MC;
    // Def. a list which contains no. of mesh points
    int m=6;
    n_MC = new int[m];
    n_MC[0] = 1000; n_MC[1] = 10000; n_MC[2] = 100000; n_MC[3] = 1000000; n_MC[4] = 10000000; n_MC[5] = 100000000;
    printf("Brute force Monte Carlo: \n");
    printf("Integral:     Variance:     Standard Deviation:      n: \n");
    for(int i=0; i<m; i++){
        start = clock();
        MonteCarloBruteForce(n_MC[i]);
        finish = clock();
        cout << "Elapsed time brute force MC: " << ((double) (finish-start)/CLOCKS_PER_SEC) << "  for n=" << n_MC[i] << endl;
    }
    printf("Improved Monte Carlo: \n");
    printf("Integral:     Variance:     Standard Deviation:      n: \n");
    for(int i=0; i<m; i++){
        start = clock();
        MonteCarlo(n_MC[i]);
        finish = clock();
        cout << "Elapsed time MC: " << ((double) (finish-start)/CLOCKS_PER_SEC) << "  for n=" << n_MC[i] << endl;
    }

    cout << "Exact answer: " << 5*pi*pi / (16*16) << endl;
    return 0;
}

void StupidIntegrationMethod_aka_Gauss_Legendre(int n){
    // Function for integrating with Gauss-Legendre
    // when using cartesian coordinates

    double a = -2, b = 2;       // Integration limits
    double pi = 4*atan(1);
    double Exact = 5*pi*pi / (16*16);
    // Using dynamic memory allocation and reserving memory for
    // vectors containing mesh point weights and function values
    double *x = new double [n];
    double *w = new double [n];

    // Calling Gauss-Legendre-function from lib.h
    // NB!: This function only calculates the 1D-case.
    gauleg(a, b, x, w, n);

    // Evaluating an integral with Gauss-Legendre is done
    // by performing a sum over a function multiplied with
    // the weights calculated by gauleg.

    double IntegralGaussLegendre = 0;
    // One for-loop for every dimension in the integral
    for(int a=0; a<n; a++){
        for(int b=0; b<n; b++){
            for(int c=0; c<n; c++){
                for(int d=0; d<n; d++){
                    for(int e=0; e<n; e++){
                        for(int f=0; f<n; f++){
                            IntegralGaussLegendre += w[a]*w[b]*w[c]*w[d]*w[e]*w[f] * Integrand(x[a],x[b],x[c],x[d],x[e],x[f]);
                        }
                    }
                }
            }
        }
    }
    printf("%f        %f        %d \n", IntegralGaussLegendre, (IntegralGaussLegendre - Exact)/Exact, n);
}
double Integrand(double x1, double y1, double z1, double x2, double y2, double z2){
    // Function for calculating integrand value at a point (x1,x2,y1,y2,z1,z2)
    // when using Gauss-Legendre and cartesian coordinates
    double alpha = 2.0;
    double ExponentialFactor1 = (-2*alpha*sqrt(x1*x1 + y1*y1 + z1*z1));
    double ExponentialFactor2 = (-2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2));
    double denominator = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
    // if-test for avoiding having denominator = 0
    if(denominator < pow(10,-6)){
        return 0;
    }
    else
        return exp(ExponentialFactor1 + ExponentialFactor2) / denominator;
}

void ABitBetterIntegrationMethod_aka_Gauss_Laguerre(int n){
    // Function for integrating with Gauss-Laguerre
    // and Gauss-Legendre when using spherical coordinates
    double pi = atan(1)*4;
    double Exact = 5*pi*pi / (16*16);
    double theta_min=0, theta_max=pi, phi_min=0, phi_max=2*pi; // Integration limits
    double alpha = 2.;

    // Distances can be integrated via Laguerre
    double *r1 = new double[n+1];
    double *wr1 = new double[n+1];

    // Angles can be integrated via Legendre
    double *theta = new double[n];
    double *wtheta = new double[n];
    double *phi = new double[n];
    double *wphi = new double[n];

    gauss_laguerre(r1, wr1, n, alpha);
    gauleg(theta_min,theta_max,theta,wtheta,n);
    gauleg(phi_min,phi_max,phi,wphi,n);

    double integral_laguerre = 0;
    for(int i=1; i<n+1; i++){
        for(int j=1; j<n+1; j++){
            for(int k=0; k<n; k++){
                for(int l=0; l<n; l++){
                    for(int m=0; m<n; m++){
                        for(int o=0; o<n; o++){
                            integral_laguerre += wr1[i]*wr1[j]*wtheta[k]*wtheta[l]*wphi[m]*wphi[o]*IntegrandSphericalCoordinates(r1[i],r1[j],theta[k],theta[l],phi[m],phi[o]);
                        }
                    }
                }
            }
        }
    }
    printf("%f        %f        %d \n", integral_laguerre, (integral_laguerre - Exact)/Exact, n);
}

double IntegrandSphericalCoordinates(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
    // Function for calculating value of integrand at (r1,r2,theta1,theta2,phi1,phi2) when using Gauss-Laguerre
    // and spherical coordinates.
    double alpha = 2;
    double ExponentialFactor = -(2*alpha - 1 )*(r1 + r2);
    double denominator = (r1*r1 + r2*r2 - 2*r1*r2*(cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*(cos(phi1 - phi2))));
    if(denominator < pow(10,-12)){
        return 0;
    }
    else
        return exp(ExponentialFactor)*sin(theta1) * sin(theta2) / sqrt(denominator); // where does the strange number come from?
}

void MonteCarloBruteForce(int n){
    double MCint = 0, MCintsqr2 = 0, temp_MCint, variance;

    long idum=-1;                       // This is the seed
    double x[6];
    double length=3.;                               // Integration limit
    double JacobiDeterminant=pow((2*length),6);     // we'll integrate over a 6D-volume box with sides equal 'length'

    // First for loop goes over every integration point
    for(int i=1; i<n; i++){
        // Second for-loop generates six random numbers to be used in integration
        for(int j=0; j<6; j++){
            // use substitution of variables, z \in [a,b] -> x \in [0,1] => z = a + (b-a)*x
            x[j] = -length+2*length*ran0(&idum);
        }
        temp_MCint = IntegrandMC(x);
        MCint += temp_MCint;
        MCintsqr2 += temp_MCint*temp_MCint;
    }

    MCint       = JacobiDeterminant*MCint/((double) n);
    MCintsqr2   = JacobiDeterminant*MCintsqr2/((double) n);
    variance    = MCintsqr2 - MCint*MCint;
    double stdev = JacobiDeterminant*sqrt(fabs(variance)/n);
    printf("%f      %f      %f                   %d \n", MCint, variance,stdev,n);
}

double IntegrandMC(double *x){
    // Calculating integrand when using brute force Monte Carlo
    double alpha = 2.0;
    double ExponentialFactor1 = (-2*alpha*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
    double ExponentialFactor2 = (-2*alpha*sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]));
    double denominator = sqrt( pow(x[0]-x[3],2) + pow(x[1]-x[4],2) + pow(x[2]-x[5],2));
    // if-test for avoiding having denominator = 0
    if(denominator < pow(10,-6)){
        return 0;
    }
    else
        return exp(ExponentialFactor1 + ExponentialFactor2) / denominator;
}

void MonteCarlo(int n){
    double pi = 4*atan(1);
    double MCint = 0, MCintsqr2 = 0, temp_MCint, variance;

    long idum=-1;                       // This is the seed
    double x[6];
    double JacobiDeterminant = 4*pow(pi,4.) *(1./16); // What about JacobiDeterminant from r1 and r2? *1.1954 ;//

    // First for loop goes over every integration point, this time with importance sampling
    for(int i=1; i<n; i++){
        // First loop is for radial part
        for(int j=0; j<2; j++){
            double y = ran0(&idum);
            x[j] = -0.25*log(1. - y); // The exponential dist. is the origin for this expression, idk why .25 is there
        }
        // Loops for the angular part, generating numbers in the interval [0,pi) and [0,2pi)
        for(int j=2; j<4; j++){
            x[j] = pi*ran0(&idum);
        }
        for(int j=4; j<6; j++){
            x[j] = 2*pi*ran0(&idum);
        }
        temp_MCint = IntegrandMC_radial(x);
        MCint += temp_MCint;
        MCintsqr2 += temp_MCint*temp_MCint;
    }

    MCint       = JacobiDeterminant*MCint/((double) n);
    MCintsqr2   = JacobiDeterminant*MCintsqr2/((double) n);
    variance    = MCintsqr2 - MCint*MCint;
    double stdev = JacobiDeterminant*sqrt(fabs(variance)/n);
    printf("%f      %f      %f                %d \n", MCint, variance,stdev,n);
}

double IntegrandMC_radial(double *x){
    // Calculating integrand when using spherical coordinates and
    // improved Monte Carlo
    double denominator = sqrt(x[0]*x[0] + x[1]*x[1] - 2*x[0]*x[1]*(cos(x[2])*cos(x[3]) + sin(x[2])*sin(x[3])*cos(x[4] - x[5]) ) );
    /*
     * x[0] = r1, x[1] = r2, x[2] = theta1, x[3] = theta2, x[4] = phi1, x[5] = phi2
     */
    return x[0]*x[0] * x[1]*x[1] * sin(x[2]) * sin(x[3]) / denominator;
    // exp(ExponentialFactor) *
}
