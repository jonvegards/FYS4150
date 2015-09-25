#include <iostream>
#include <armadillo>
#include <string>
#include "lib.h"
using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void rotateAmatrix(double &, double &, int &, int &, int , mat &, mat &R);
void findSinCos(mat &, int &, int &, double &, double &);
void printMatlabMatrix(string name, mat &A);
void JacobiRotation(int, double rho_max, mat A);

int main()
{
    cout.precision(16); // Setting how many digits the program should print
    //======================================================
    // Defining n and calling on the rotation-function
    //======================================================
    int n = 200;
    double rho_max = 5.2;

    //======================================================
    // Defining parameters, step length, ...
    //======================================================
    double e, h, h_temp;

    h = rho_max / (n+2);
    e = -1 / (h*h);
    h_temp = 2 / (h*h);

    //======================================================
    // Making the rhoes and potential vector
    //======================================================
    vec rho(n+2);
    for(int i=0; i<n+2; i++){
        rho(i) += h*i;
    }

    vec v = rho % rho; // % is the operator for inner product

    //======================================================
    // Declaring matrix A w/elements (see eq. (2) in exercises)
    //======================================================
    mat A(n,n);
    for(int i=0; i<n; i++){
        A(i,i) = h_temp + v(i+1);
    }
    A.diag(1) += e;
    A.diag(-1) += e;

    mat A_armadillo = A;

    JacobiRotation(n, rho_max, A);

    //======================================================
    // Using Armadillo to solve the eigenvalue equation
    //======================================================
    vec eigval = eig_sym( A_armadillo );

    //======================================================
    // Printing out results
    //======================================================
    cout << "Calculation by Armadillo:" << endl;
    cout << "lambda0: " << eigval(0) << "  lambda1: " << eigval(1) << "  lambda2: " << eigval(2) << endl;
}

void JacobiRotation(int n, double rho_max, mat A){

    //======================================================
    // Finding biggest element in A
    //======================================================
    int k,l;
    double max_A = 0;
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);

    //======================================================
    // Setting limits on the while-loop's duration and
    // accuracy
    //======================================================
    int iterations=0;
    int max_number_iterations = pow(10,5);
    double epsilon = pow(10,-8);

    //======================================================
    // Printing A for checking with MATLAB
    //======================================================
    //printMatlabMatrix("A", A);

    //======================================================
    // Calculating the equation
    //======================================================
    mat R = eye<mat>(n,n);
    while( fabs(max_A) > epsilon && (double) iterations < max_number_iterations){
        double c=0, s=0;
        findSinCos(A, k, l, s, c);
        rotateAmatrix(s,c,k,l,n,A,R);

        // Find max off-diagonal matrix element
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        //cout << "max_A: " << max_A << endl;
        iterations++;
    }

    //======================================================
    // Preparing results for printing
    //======================================================
    //R.print();
    vec diagonalA(n);
    diagonalA = sort(A.diag()); // sort(A) reorder the elements so they lowest values come first

    //======================================================
    // Printing out results
    //======================================================
    cout << "n: " << n << endl;
    cout << "rho_max: " << rho_max << endl;
    cout << "# of iterations: " << iterations << endl;
    cout << "lambda0: " << diagonalA(0) << "  lambda1: " << diagonalA(1) << "  lambda2: " << diagonalA(2) << endl;
}

void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n){
    // Find max off-diagonal matrix element
    max_A = 0.0;
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if( fabs(A(i,j)) > max_A){
                k = i;
                l = j;
                max_A = fabs(A(i,j));
            }
        }
    }
}

void findSinCos(mat &A, int &k, int &l, double &s, double &c){
    double tau, t, temp_sqrt, t1, t2;

    tau = -( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
    temp_sqrt = sqrt( 1 + tau*tau );
    t1 = tau + temp_sqrt;
    t2 = tau - temp_sqrt;

    if(fabs(t1) < fabs(t2)){
        t = t1;
    }
    else{
        t = t2;
    }

    c = 1 / sqrt( 1 + t*t);
    s = c*t;
}

void rotateAmatrix(double &s, double &c, int &k, int &l, int n, mat &A, mat &R){
    // Rotate matrix A
    for(int i=0; i<n; i++){
        if(i != k && i != l){
            A(i,i) = A(i,i);
            double a_il = A(i,l);
            double a_ik = A(i,k);
            A(i,k) = a_ik*c - a_il*s;
            A(i,l) = a_il*c + a_ik*s;
            // A is symmetric so we can only copy the already calculated values
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }
        // Finding new eigenvectors
        double r_ik = R(i,k);
        double r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }

    // Rotating the diagonal elements of A
    double a_kk = A(k,k);
    double a_ll = A(l,l);
    A(k,k) = a_kk*c*c - 2*A(k,l)*s*c + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*A(k,l)*s*c + a_kk*s*s;
    A(k,l) = 0; // Why is this zero?
    A(l,k) = 0;
}

//======================================================
// Function for printing matrix A so it can be copy-
// pasted into MATLAB for testing of this program's
// capability to find the eigenvalues.
//======================================================
void printMatlabMatrix(string name, mat &A) {
    cout << name << " = [";
    for(int i=0; i<A.n_rows; i++) {
        for(int j=0; j<A.n_cols; j++) {
            cout << A(i,j) << " ";
        }
        cout << "; ";
    }
    cout << "];" << endl;
    cout << "[U,V] = eig("<<name<<")" << endl;
}
