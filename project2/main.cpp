#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void rotateAmatrix(double &, double &, int &, int &, int , mat &, mat &);
void findSinCos(mat &, int &, int &, double &, double &);
void printMatlabMatrix(string name, mat &);
void JacobiRotation(int, mat & , uvec &, mat &);
void TwoElectronCase();
void TautAnalyticalSolution(int , double , vec &);
void SavingResultsToFile(string , string , string , vec , vec );

// Functions for testing the program
void RunAllTests();
void TwoByTwoMatrixTest();
void MaxElementTest();
void EigenvalueTest();

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
    int n = 250;
    double rho_max = 12.0; // The smallest omega-values required rho_max=18 to get a stable solution

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

//    cout << "Single Electron Case" << endl;
//    JacobiRotation(n, A_SEC, diagonalA_SEC, R_SEC);

//    cout << "n: " << n << endl;
//    cout << "rho_max: " << rho_max << endl;
//    cout << "lambda0: " << A.diag()(diagonalA(0)) << "  lambda1: " << A.diag()(diagonalA(1)) << "  lambda2: " << A.diag()(diagonalA(2)) << endl;

    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    //------------------------------------------------------
    mat eigvec1;
    vec eigval1;

//    eig_sym(eigval1, eigvec1, A_SEC_armadillo);

//    cout << "Calculation by Armadillo:" << endl;
//    cout << "lambda0: " << eigval1(0) << "  lambda1: " << eigval1(1) << "  lambda2: " << eigval1(2) << endl;

    //============================================================================================================
    // Two electron case (TEC)
    //============================================================================================================
    double omega_r = .25; // set this to 0.25 for checking with analytical sol.
    vec TwoElectronPotential = omega_r*omega_r*(rho % rho) + 1 / rho;
    mat A_TEC(n,n);
    for(int i=0; i<n; i++){
        A_TEC(i,i) = h_temp + TwoElectronPotential(i+1);
    }
    A_TEC.diag(1) += e;
    A_TEC.diag(-1) += e;
    mat A_TEC_armadillo = A_TEC;
    uvec diagonalA_TEC;
    mat R_TEC = eye<mat>(n,n);

//    cout << " " << endl;
    cout << "Two electron case:" << endl;
    JacobiRotation(n, A_TEC, diagonalA_TEC, R_TEC);

    cout << "n: " << n << endl;
    cout << "rho_max: " << rho_max << endl;
    cout << "lambda0: " << A_TEC.diag()(diagonalA_TEC(0)) << "  lambda1: " << A_TEC.diag()(diagonalA_TEC(1)) << "  lambda2: " << A_TEC.diag()(diagonalA_TEC(2)) << endl;


    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    // with Coulomb-repulsion
    //------------------------------------------------------
//    mat eigvec2;
//    vec eigval2;

//    eig_sym(eigval2, eigvec2, A_TEC_armadillo);

//    cout << "Calculation by Armadillo:" << endl;
//    cout << "lambda0: " << eigval2(0) << "  lambda1: " << eigval2(1) << "  lambda2: " << eigval2(2) << endl;

    //============================================================================================================
    // Analytical solution
    //============================================================================================================
    vec Psi_anal;
    TautAnalyticalSolution(n, h, Psi_anal);
    //Psi_anal.print();

    //============================================================================================================
    // Saving results to MATLAB-files
    //============================================================================================================
    // First our result from Jacobi's method TEC and SEC
//    SavingResultsToFile("psi_and_rho_SEC.m", "R", "rho", R_SEC.col(diagonalA_SEC(0)), rho);

    SavingResultsToFile("psi_and_rho_TEC1.m", "R", "rho", R_TEC.col(diagonalA_TEC(0)), rho);

    // R.col(diagonalA(0)) is the eigenvector belonging to
    // the lowest eigenvalue
    // Remember to not plot first and last element in rho-vector (boundary conditions)

    // Saving analytical solution
    SavingResultsToFile("psi_and_rho_anal.m", "psi", "rho", Psi_anal, rho);

    // Armadillo's solutions TEC and SEC
    //SavingResultsToFile("psi_and_rho_armadillo_SEC.m", "psi", "rho", eigvec1.col(0), rho);
    //SavingResultsToFile("psi_and_rho_armadillo_TEC.m", "psi", "rho", eigvec2.col(0), rho);

    //============================================================================================================
    // Running unit tests
    //============================================================================================================
    //RunAllTests();
}

void JacobiRotation(int n, mat &A, uvec &diagonalA, mat &R){

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
    cout << "# of iterations: " << iterations << endl;
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

//======================================================
// Unit tests
//======================================================

void RunAllTests(){
    TwoByTwoMatrixTest();
    MaxElementTest();
    EigenvalueTest();
}

void TwoByTwoMatrixTest(){
    int n = 2;
    mat A(n,n);
    A(0,0) = 1;
    A(1,1) = 2;
    A.diag(1) += 1;
    A.diag(-1) += 1;

    mat R = eye<mat>(n,n);
    uvec diagonalA;

    cout << "Testing function JacobiRotation with a (2x2)-matrix A =" << endl;
    A.print();
    cout << "After running the function A =" << endl;
    JacobiRotation(n, A, diagonalA, R);
    cout << "A = " << endl;
    A.print();
}

void MaxElementTest(){
    int k,l, n = 3;
    double max_A;
    mat A = zeros<mat>(n,n);
    A.diag() += 1;
    A(0,2) = 1;
    A(0,2) = 3;
    A(1,2) = -100;
    A(1,0) = A(2,0) = A(2,1) = 10000; // findMax...-function searches only in upper triangular
    cout << "Test of function findMaximumElementOnNonDiagonal. With the matrix"<< endl;
    cout << "A = " << endl;
    A.print();
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
    cout << "Function found max element to be: " << A(k,l) << endl;
    cout << "The function searches only the upper triangular part, so it should give us max(A) = " << A(1,2) << endl;

//    mat B(n,n);
//    B(0,0) = 1;
//    B(1,1) = 2;
//    B.diag(1) += 1000;
//    B.diag(-1) += 0.2;
//    findMaximumElementOnNonDiagonal(B, o, p, max_A, n)
}

void EigenvalueTest(){
    // Function for testing if the program calculates the correct eigenvalues
    // for quantum harmonic oscillator. Analytical eigenvalues are
    // lambda0 = 3, lambda1 = 7, lambda2 = 11.

    int n = 200;
    double rho_max = 5.0;

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

    // The potential is rho^2
    vec v = rho % rho; // % is the operator for inner product
    mat A(n,n);
    for(int i=0; i<n; i++){
        A(i,i) = h_temp + v(i+1);
    }
    A.diag(1) += e;
    A.diag(-1) += e;

    uvec diagonalA; // vector to contain eigenvalue-indices
    mat R = eye<mat>(n,n); // matrix to contain eigenvectors
    cout << "Checking if the program gives the correct eigenvalues" << endl;
    cout << "Analytical solutions: lambda0 = 3, lambda1 = 7, lambda2 = 11" << endl;
    cout << "The calculated solutions are: " << endl;
    JacobiRotation(n, A, diagonalA, R);
    cout << "# of mesh points n: " << n << endl;
    cout << "lambda0: " << A.diag()(diagonalA(0)) << "  lambda1: " << A.diag()(diagonalA(1)) << "  lambda2: " << A.diag()(diagonalA(2)) << endl;
}
