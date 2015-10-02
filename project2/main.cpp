#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void rotateAmatrix(double &, double &, int &, int &, int , mat &, mat &);
void findSinCos(mat &, int &, int &, double &, double &);
void printMatlabMatrix(string name, mat &);
void JacobiRotation(int, double , mat , uvec &, mat &);
void TwoElectronCase();
void TautAnalyticalSolution(int , double , vec &);
void SavingResultsToFile(string , string , string , vec , vec );

//======================================================
// This program calculates a matrix eigenvalue equation
// by using Jacobi's method. It should work for every
// (n x n)-matrix. The eigenvectors together with the
// distance vector are written to files (almost) ready to
// be run i MATLAB. The cases we're looking at here is
// the Schrodinger eq. for a harmonic oscillator with one
// and two electrons, the only difference is the potential
// used. The analytical solution is also calculated.
//======================================================

int main()
{
    cout.precision(6); // Setting how many digits the program should print
    //======================================================
    // Defining n and calling on the rotation-function
    int n = 500;
    double rho_max = 18.0; // The smallest omega-values required big rho to get a stable solution

    //======================================================
    // Defining parameters, step length, ...
    double e, h, h_temp;
    h = rho_max / (n+2);
    e = -1 / (h*h);
    h_temp = 2 / (h*h);

    //======================================================
    // Making the rho-vector
    vec rho(n+2);
    for(int i=0; i<n+2; i++){
        rho(i) += h*i;
    }
    rho(0) = 0;

    //============================================================================================================
    // Single electron case (SEC)
    //============================================================================================================
    // The potential is rho^2
    vec v = rho % rho; // % is the operator for inner product
    mat A_SEC(n,n);
    for(int i=0; i<n; i++){
        A_SEC(i,i) = h_temp + v(i+1);
    }
    A_SEC.diag(1) += e;
    A_SEC.diag(-1) += e;

    mat A_SEC_armadillo = A_SEC; // Save matrix for using with Armadillo
    uvec diagonalA_SEC;
    mat R_SEC = eye<mat>(n,n);

    cout << "Single Electron Case" << endl;
    //JacobiRotation(n, rho_max, A_SEC, diagonalA_SEC, R_SEC);

    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    //------------------------------------------------------
    mat eigvec1;
    vec eigval1;
    //eig_sym(eigval1, eigvec1, A_SEC_armadillo);

    // Printing out results from Armadillo-solver
    cout << "Calculation by Armadillo:" << endl;
    //cout << "lambda0: " << eigval1(0) << "  lambda1: " << eigval1(1) << "  lambda2: " << eigval1(2) << endl;

    //============================================================================================================
    // Two electron case (TEC)
    //============================================================================================================
    double omega_r = 5;
    vec TwoElectronPotential = omega_r*(rho % rho) + 1 / rho;
    mat A_TEC(n,n);
    for(int i=0; i<n; i++){
        A_TEC(i,i) = h_temp + TwoElectronPotential(i+1);
    }
    A_TEC.diag(1) += e;
    A_TEC.diag(-1) += e;
    mat A_TEC_armadillo = A_TEC;
    uvec diagonalA_TEC;
    mat R_TEC = eye<mat>(n,n);

    cout << " " << endl;
    cout << "Two electron case:" << endl;
    //JacobiRotation(n, rho_max, A_TEC, diagonalA_TEC, R_TEC);


    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    // with Coulomb-repulsion
    //------------------------------------------------------
    mat eigvec2;
    vec eigval2;
    eig_sym(eigval2, eigvec2, A_TEC_armadillo);

    // Printing out results from Armadillo-solver
    cout << "Calculation by Armadillo:" << endl;
    cout << "lambda0: " << eigval2(0) << "  lambda1: " << eigval2(1) << "  lambda2: " << eigval2(2) << endl;

    //============================================================================================================
    // Analytical solution
    //============================================================================================================
    vec Psi_anal;
    //TautAnalyticalSolution(n, h, Psi_anal);
    //Psi_anal.print();

    //============================================================================================================
    // Saving results to file
    //============================================================================================================
    // First our result from Jacobi's method (TEC)
//    SavingResultsToFile("psi_and_rho_SEC1.m", "R", "rho", R_SEC.col(diagonalA_SEC(0)), rho);
    // First our result from Jacobi's method (TEC)
//    SavingResultsToFile("psi_and_rho_TEC1.m", "R", "rho", R_TEC.col(diagonalA_TEC(0)), rho);

    // R.col(diagonalA(0)) is the eigenvector belonging to
    // the lowest eigenvalue
    // Remember to not plot first and last element in rho-vector (boundary conditions)

    // Analytical solution
    //SavingResultsToFile("psi_and_rho_anal.m", "psi", "rho", Psi_anal, rho);

    // Armadillo's solution SEC
    //SavingResultsToFile("psi_and_rho_armadillo_SEC.m", "psi", "rho", eigvec1.col(0), rho);
    // Armadillo's solution TEC
    SavingResultsToFile("psi_and_rho_armadillo_TEC1.m", "psi", "rho", eigvec2.col(0), rho);
}

void JacobiRotation(int n, double rho_max, mat A, uvec &diagonalA, mat &R){

    //======================================================
    // Finding biggest element in A
    //======================================================
    int k,l;
    double max_A = 0;
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);

    //======================================================
    // Setting limits on the while-loop's duration and
    // accuracy
    int iterations=0;
    int max_number_iterations = pow(10,6);
    double epsilon = pow(10,-8);

    //======================================================
    // Printing A for checking with MATLAB
    //printMatlabMatrix("A", A);

    //======================================================
    // Calculating the equation
    //mat R = eye<mat>(n,n);

    while( fabs(max_A) > epsilon && (double) iterations < max_number_iterations){
        double c=0, s=0;
        findSinCos(A, k, l, s, c);
        rotateAmatrix(s,c,k,l,n,A,R);
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        iterations++;
    }

    //======================================================
    // Preparing results for printing
    diagonalA = sort_index(A.diag());

    //======================================================
    // Printing out results
    cout << "n: " << n << endl;
    cout << "rho_max: " << rho_max << endl;
    cout << "# of iterations: " << iterations << endl;
    cout << "lambda0: " << A.diag()(diagonalA(0)) << "  lambda1: " << A.diag()(diagonalA(1)) << "  lambda2: " << A.diag()(diagonalA(2)) << endl;

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
    double tau, t, temp_sqrt;
    tau = -( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
    temp_sqrt = sqrt( 1 + tau*tau );
    if(tau < 0){
        t = tau + temp_sqrt;
    }
    else{
        t = tau - temp_sqrt;
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

void printMatlabMatrix(string name, mat &A) {
    //======================================================
    // Function for printing (m x n)-matrix A so it can be copy-
    // pasted into MATLAB for testing of this program's
    // capability to find the eigenvalues.
    //======================================================
    cout << name << " = [";
    for(int i=0; i<int(A.n_rows); i++) {
        for(int j=0; j<int(A.n_cols); j++) {
            cout << A(i,j) << " ";
        }
        cout << "; ";
    }
    cout << "];" << endl;
    cout << "[U,V] = eig("<<name<<")" << endl;
}

void TautAnalyticalSolution(int n, double h, vec &Psi_anal){
    //======================================================
    // Calculating Taut's analytical solution
    //======================================================
    vec r(n+2);
    vec phi(n+2);
    int l=0;

    for(int i=0; i<n+2; i++){
        r(i) += i*h;
    }

    r(0) = 0;
    phi = pow(r, 1) % exp(- pow(r, 2) / (8*(l+1)) ) % ( 1 + r / (2*(l+1) ) );
    //phi = exp(- pow(r, 2));
    //phi(0) = r(0) = 0;

    vec ksi = exp(-r%r);
    Psi_anal = phi; //% ksi;
}

void SavingResultsToFile(string FileName, string namevec1, string namevec2, vec vec1, vec vec2){
    ofstream myfile;
    myfile.open (FileName);
    myfile << namevec1 << "= [";
    for (int i=0; i<int(vec1.n_rows); i++){
        myfile << vec1(i) << ", ";
    }
    myfile << "];" << endl;

    myfile << namevec2 << "= [";
    for (int i=0; i<int(vec2.n_rows); i++){
        myfile << vec2(i) << ", ";
    }
    myfile << "];" << endl;
    myfile << "plot("<< namevec2 <<","<< namevec1 << ".*"<< namevec1 << ")" << endl;
    myfile << "n=" << vec2.n_rows << endl;
    myfile.close();
}
