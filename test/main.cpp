#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;


int main(){
  mat A = randu<mat>(5,5);
  vec b = randu<vec>(5);
  A.print("A =");
  b.print("b=");


  // some simple matrix operations
  // determinant
  cout << "det(A) = " << det(A) << endl;
  // inverse
  cout << "inv(A) = " << endl << inv(A) << endl;


  // solve Ax = b
  vec x = solve(A,b);
  // print x
  x.print("x=");
  // find LU decomp of A, P is the permutation matrix
  mat L, U, P;
  lu(L,U,P, A);
  // print l
  L.print("L=");
  // print U
  U.print("U=");
    return 0;
  }
