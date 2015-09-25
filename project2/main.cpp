#include <iostream>
#include <armadillo>
#include <string>
using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void rotateAmatrix(double &, double &, int &, int &, int , mat &);
void findSinCos(mat &, int &, int &, double &, double &);
void printMatlabMatrix(string name, mat &A);

int main()
{
    double e, h, h_temp, rho_max, epsilon;
    int n=3;
    // n is here the dimension of the matrix we're going to use
    // no. of steps is then (n+2)

    mat A(n,n);
    vec rho(n+2);
    vec v(n+2);

    // Defining parameters
    epsilon = pow(10,-8);

    rho_max = 8.;
    h = rho_max / (n+2);
    e = -1 / (h*h);
    h_temp = 2 / (h*h);

    // Making the rhoes
    for(int i=0; i<n+2; i++){
        rho(i) += h*i;
    }

    // Vector w/discretized potential
    // % is the operator for inner product
    v = rho % rho;

    // Filling matrix A w/elements (see eq. (2) in exercises)
    for(int i=0; i<n; i++){
        A(i,i) = h_temp + v(i+1);
    }
    A.diag(1) += e;
    A.diag(-1) += e;

    // Starting w/A(1,2)
    int k = 0;
    int l = 1;
    double max_A = 0;
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);

    int iterations=0;
    int max_number_iterations = 1000;
    A.print();

    while( fabs(max_A) > epsilon && (double) iterations < max_number_iterations){
        double c=0, s=0;
        findSinCos(A, k, l, s, c);
        rotateAmatrix(s,c,k,l,n,A);

        // Find max off-diagonal matrix element
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        cout << "max_A: " << max_A << endl;
        iterations++;
    }
    cout << "# of iterations: " << iterations << endl;
    // A.print();
    printMatlabMatrix("A", A);
}


void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n){
    // Find max off-diagonal matrix element
    max_A = 0.0;
    for(int i=0; i<n; i++){
        //cout << max_A << endl;
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

void rotateAmatrix(double &s, double &c, int &k, int &l, int n, mat &A){
    for(int i=0; i<n; i++){
        if(i != k && i != l){
            A(i,i) = A(i,i);
            double a_il = A(i,l);
            double a_ik = A(i,k);
            A(i,k) = a_ik*c - a_il*s;
            A(i,l) = a_il*c + a_ik*s;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }
    }

    A(k,k) = A(k,k)*c*c - 2*A(k,l)*s*c + A(l,l)*s*s;
    A(l,l) = A(l,l)*c*c + 2*A(k,l)*s*c + A(k,k)*s*s;
    A(k,l) = 0;
    A(l,k) = 0;
}

void printMatlabMatrix(string name, mat &A) {
    cout << name << " = [";
    for(int i=0; i<A.n_rows; i++) {
        for(int j=0; j<A.n_cols; j++) {
            cout << A(i,j) << " ";
        }
        cout << "; ";
    }
    cout << "];" << endl;
    cout << "eig("<<name<<")" << endl;
}
