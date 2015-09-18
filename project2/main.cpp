#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &, int &, int &, double &, int);
void maxOfT1andT2(double &, double &, double &);

int main()
{
    int i, n=10, k, l;
    double e, h, h_temp, rho_max, epsilon, tau, t1, t2, c, t, s;

    mat A(n,n), B(n,n), S(n,n);
    vec rho(n);
    vec v(n);

    // Defining parameters
    epsilon = pow(10,-8);
    cout << epsilon << endl;
    rho_max = 8.;
    h = rho_max / n;
    e = -1 / (h*h);
    h_temp = 2 / (h*h);

    // Making the rhoes
    for(i=0; i<n; i++){
        rho(i) += h*i;
    }

    // Vector w/discretized potential
    v = rho % rho;
    v.print();

    // Filling matrix A w/elements (see eq. (2) in exercises)
    for(i=0; i<n; i++){
        A(i,i) = h_temp + v(i);
    }
    A.diag(1) += e;
    A.diag(-1) += e;

    // Starting w/A(1,2)
    k = 0;
    l = 1;

    double max_A = A(k,l);
    double temp_sqrt;

    while( fabs(max_A) > epsilon){
        tau = -( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
        temp_sqrt = sqrt( 1 + tau*tau );
        t1 = tau + temp_sqrt;
        t2 = tau - temp_sqrt;
        maxOfT1andT2(t1,t2,t); // max of t1 vs. t2
        c = 1 / sqrt( 1 + t*t);
        s = c*t;
        S.diag() += 1.;
        S(k,k) = c;
        S(l,l) = c;
        S(k,k) = s;
        S(l,l) = -s;
        A = S.t() * A * S;
        // Find max off-diagonal matrix element
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        //cout << max_A << endl;
    }

    //A.print();

}
void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n){
    // Find max off-diagonal matrix element
    int i, j;
    for(i=1; i<n; i++){
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

void maxOfT1andT2(double &t1, double &t2, double &t){
    if(fabs(t1) < fabs(t2)){
        t = t1;
    }
    else{
        t = t2;
    }
}

