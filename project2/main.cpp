#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void rotateAmatrix(double &, double &, int &, int &, int , mat &, mat &);
void findSinCos(mat &, int &, int &, double &, double &);
void printMatlabMatrix(string name, mat &);
void JacobiRotation(int, mat & , uvec &, mat &, int &);
void TwoElectronCase();
void TautAnalyticalSolution(int , double , vec &);
void SavingResultsToFile(string , string , string , vec , vec );
void SavingParametersToFile(int, double, mat, uvec, mat, uvec, vec, vec, int, int, double);
void MainCall();

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
    //====================================================
    // Calling on function to do the calculations,
    // this must be commented out if you want to run
    // RunAllTests();
    MainCall();

    //====================================================
    // Running unit tests
    //RunAllTests();

    // Due to some memory problems we cannot run MainCall();
    // and RunAllTests(); simultaneously. Why?
}

void MainCall(){
    cout.precision(6); // Setting how many digits the program should print out
    //------------------------------------------------------
    int n = 100;
    double rho_max = 12.0; // The smallest omega-values required rho_max=18 to get a stable solution
    // NB!: When using n=2 you must not print out the three lowest eigenvalues
    // since it's only two eigenvalues.

    //------------------------------------------------------
    double h = rho_max / (n+2); // Step length
    double e = -1 / (h*h);      // Elements on upper/lower diagonal
    double h_temp = 2 / (h*h);  // Temporary variable to save FLOPS
    // Note that we use (n+2) steps in rho, this is because
    // we have to omit the end points (since they are zero)
    // when plugging the initial conditions into the matrix.

    //------------------------------------------------------
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
    mat A_SEC = zeros<mat>(n,n);
    // Filling matrix A_SEC on the diagonal with the potential terms and terms from three point-formula
    for(int i=0; i<n; i++){
        A_SEC(i,i) = h_temp + v(i+1);
    }
    A_SEC.diag(1) += e; // Filling the diaognals over and under the main diagonal
    A_SEC.diag(-1) += e;

    mat A_SEC_armadillo = A_SEC; // Save matrix for using with Armadillo
    uvec diagonalA_SEC;
    mat R_SEC = eye<mat>(n,n);
    int iterations_SEC;
    JacobiRotation(n, A_SEC, diagonalA_SEC, R_SEC, iterations_SEC);

    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    //------------------------------------------------------
    mat eigvec1;
    vec eigval1;

    eig_sym(eigval1, eigvec1, A_SEC_armadillo);

    //============================================================================================================
    // Two electron case (TEC)
    //============================================================================================================

    double omega_r = 1; // set this to 0.25 for checking with analytical sol.
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
    int iterations_TEC;
    JacobiRotation(n, A_TEC, diagonalA_TEC, R_TEC, iterations_TEC);



    //------------------------------------------------------
    // Using Armadillo to solve the eigenvalue equation
    // with Coulomb-repulsion
    //------------------------------------------------------

    mat eigvec2;
    vec eigval2;

    eig_sym(eigval2, eigvec2, A_TEC_armadillo);

    //============================================================================================================
    // Analytical solution
    //============================================================================================================
    vec Psi_anal;
    TautAnalyticalSolution(n, h, Psi_anal);

    //============================================================================================================
    // Saving results to MATLAB-files
    //============================================================================================================
    // First our result from Jacobi's method TEC and SEC
    SavingResultsToFile("psi_and_rho_SEC.m", "R", "rho", R_SEC.col(diagonalA_SEC(0)), rho);
    SavingResultsToFile("psi_and_rho_TEC.m", "R", "rho", R_TEC.col(diagonalA_TEC(0)), rho);

    // R.col(diagonalA(0)) is the eigenvector belonging to
    // the lowest eigenvalue
    // Remember to not plot first and last element in rho-vector (due to boundary conditions)

    // Saving analytical solution
    SavingResultsToFile("psi_and_rho_anal.m", "psi", "rho", Psi_anal, rho);

    // Armadillo's solutions TEC and SEC
    SavingResultsToFile("psi_and_rho_armadillo_SEC.m", "psi", "rho", eigvec1.col(0), rho);
    SavingResultsToFile("psi_and_rho_armadillo_TEC.m", "psi", "rho", eigvec2.col(0), rho);

    //============================================================================================================
    // Saving used parameters to file
    //============================================================================================================
    SavingParametersToFile(n, rho_max, A_TEC, diagonalA_TEC, A_SEC, diagonalA_SEC, eigval1, eigval2, iterations_SEC, iterations_TEC, omega_r);

}

void JacobiRotation(int n, mat &A, uvec &diagonalA, mat &R, int &iterations){

    //------------------------------------------------------
    // Finding biggest element in A
    //------------------------------------------------------
    int k,l;
    double max_A = 0;
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);

    //------------------------------------------------------
    // Setting limits on the while-loop's duration and
    // accuracy
    iterations=0;
    int max_number_iterations = pow(10,6);
    double epsilon = pow(10,-9);

    //------------------------------------------------------
    // Printing A for checking with MATLAB
    //printMatlabMatrix("A", A);

    //------------------------------------------------------
    // Calculating the equation
    
    while( fabs(max_A) > epsilon && (double) iterations < max_number_iterations){
        double c=0, s=0;
        findSinCos(A, k, l, s, c);
        rotateAmatrix(s,c,k,l,n,A,R);
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        iterations++;

    }

    //------------------------------------------------------
    // Preparing results for printing
    diagonalA = sort_index(A.diag());
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
    Psi_anal = pow(r, 1) % exp(- pow(r, 2) / (8*(l+1)) ) % ( 1 + r / (2*(l+1) ) );
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

void SavingParametersToFile(int n, double rho_max, mat A_TEC, uvec diagonalA_TEC, mat A_SEC, uvec diagonalA_SEC, vec eigval1, vec eigval2, int iterations_SEC, int iterations_TEC, double omega_r){
    ofstream myfile;
    myfile.open ("MainProgramRunParameters.txt");
    myfile << "Single Electron Case" << endl;
    myfile << "n: " << n << endl;
    myfile << "# of iterations: " << iterations_TEC << endl;
    myfile << "rho_max: " << rho_max << endl;
    myfile << "lambda0: " << A_SEC.diag()(diagonalA_SEC(0)) << "  lambda1: " << A_SEC.diag()(diagonalA_SEC(1)) << "  lambda2: " << A_SEC.diag()(diagonalA_SEC(2)) << endl;
    myfile << "Calculation by Armadillo:" << endl;
    myfile << "lambda0: " << eigval1(0) << "  lambda1: " << eigval1(1) << "  lambda2: " << eigval1(2) << endl;
    myfile << " " << endl;
    myfile << "Two electron case:" << endl;
    myfile << "n: " << n << endl;
    myfile << "# of iterations: " << iterations_SEC << endl;
    myfile << "rho_max: " << rho_max << endl;
    myfile << "omega_r: " << omega_r << endl;
    myfile << "lambda0: " << A_TEC.diag()(diagonalA_TEC(0)) << "  lambda1: " << A_TEC.diag()(diagonalA_TEC(1)) << "  lambda2: " << A_TEC.diag()(diagonalA_TEC(2)) << endl;
    myfile << "Calculation by Armadillo:" << endl;
    myfile << "lambda0: " << eigval2(0) << "  lambda1: " << eigval2(1) << "  lambda2: " << eigval2(2) << endl;
    myfile << "Difference between Jacobi's method and Armadillo: " << eigval2(0) - A_TEC.diag()(diagonalA_TEC(0)) << endl;
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
    int iterations22;

    cout << "Testing function JacobiRotation with a (2x2)-matrix A =" << endl;
    A.print();
    cout << "After running the function A =" << endl;
    JacobiRotation(n, A, diagonalA, R, iterations22);
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
}

void EigenvalueTest(){
    // Function for testing if the program calculates the correct eigenvalues
    // for quantum harmonic oscillator. Analytical eigenvalues are
    // lambda0 = 3, lambda1 = 7, lambda2 = 11.

    int nTest = 200;
    double rho_maxTest = 5.0;

    //------------------------------------------------------
    // Defining parameters, step length, ...

    double hTest = rho_maxTest / (nTest+2);
    double eTest = -1 / (hTest*hTest);
    double h_tempTest = 2 / (hTest*hTest);

    //------------------------------------------------------
    // Making the rho-vector
    vec rhoTest(nTest+2);
    for(int i=0; i<nTest+2; i++){
        rhoTest(i) += hTest*i;
    }
    rhoTest(0) = 0;

    // The potential is rho^2
    vec vTest = rhoTest % rhoTest;
    mat ATest(nTest,nTest);
    for(int i=0; i<nTest; i++){
        ATest(i,i) = h_tempTest + vTest(i+1);
    }
    ATest.diag(1) += eTest;
    ATest.diag(-1) += eTest;
    int iterationsTest;
    uvec diagonalATest; // vector to contain eigenvalue-indices
    mat RTest = eye<mat>(nTest,nTest); // matrix to contain eigenvectors
    cout << "Checking if the program gives the correct eigenvalues" << endl;
    cout << "Analytical solutions: lambda0 = 3, lambda1 = 7, lambda2 = 11" << endl;
    cout << "The calculated solutions are: " << endl;
    JacobiRotation(nTest, ATest, diagonalATest, RTest, iterationsTest);
    cout << "# of mesh points n: " << nTest << endl;
    cout << "# of iterations: " << iterationsTest << endl;
    cout << diagonalATest(2) << endl;
    cout << "lambda0: " << ATest.diag()(diagonalATest(0)) << "  lambda1: " << ATest.diag()(diagonalATest(1)) << "  lambda2: " << ATest.diag()(diagonalATest(2)) << endl;
}
