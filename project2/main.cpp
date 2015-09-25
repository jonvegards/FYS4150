#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void minOfT1andT2(double &, double &, double &, double &, double &);

int main()
{
    int i, k, l;
    double e, h, h_temp, rho_max, epsilon, tau, t1, t2, c, t, s, z=0;
    int n=2;
    // n is here the dimension of the matrix we're going to use
    // no. of steps is then (n+2)

    mat A(n,n), B(n,n), S(n,n);
    vec rho(n+2);
    vec v(n+2);

    // Defining parameters
    epsilon = pow(10,-8);

    rho_max = 8.;
    h = rho_max / (n+2);
    e = -1 / (h*h);
    h_temp = 2 / (h*h);

    // Making the rhoes
    for(i=0; i<n+2; i++){
        rho(i) += h*i;
    }

    // Vector w/discretized potential
    // % is the operator for inner product
    v = rho % rho;
    v.print();
    rho.print();

    // Filling matrix A w/elements (see eq. (2) in exercises)
    for(i=0; i<n; i++){
        A(i,i) = h_temp + v(i+1);
    }
    A.diag(1) += e;
    A.diag(-1) += e;

    // Starting w/A(1,2)
    k = 0;
    l = 1;

    double max_A = A(k,l);
    double temp_sqrt;

    /*
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
    cout << max_A << endl;
    A.print();*/

    while( fabs(max_A*max_A) > epsilon){
        tau = -( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
        temp_sqrt = sqrt( 1 + tau*tau );
        t1 = tau + temp_sqrt; // There should be a minus sign here, right?
        t2 = tau - temp_sqrt; // But then it doesn't workkkkk :---(
        minOfT1andT2(t1,t2,t, c, s); // min of t1 vs. t2
        cout << t <<"  " << t1 <<"  " << t2 << endl;
        S.diag() += 1.;
        S(k,k) = c;
        S(l,l) = c;
        S(k,l) = s;
        S(l,k) = -s;
        S.print();
        A = S.t() * A * S;
        // Find max off-diagonal matrix element
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        cout << max_A << endl;
        z += 1;
    }
    A.print();
    cout << z << endl;
}


void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n){
    // Find max off-diagonal matrix element
    int i=0, j=0;
    max_A = 0.0;
    for(i=0; i<n; i++){
        //cout << max_A << endl;
        for(j=i+1; j<n; j++){
            if( fabs(A(i,j)) > max_A){
                k = i;
                l = j;
                max_A = A(i,j);
            }
        }
    }
}

void minOfT1andT2(double &t1, double &t2, double &t, double &c, double &s){
    if(fabs(t1) < fabs(t2)){
        t = t1;
    }
    else{
        t = t2;
    }
    c = 1 / sqrt( 1 + t*t);
    s = c*t;
}

