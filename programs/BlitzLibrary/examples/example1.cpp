//     Simple test case of matrix operations
//     using Blitz++
#include <blitz/array.h> 
#include <iostream>
using namespace std;
using namespace blitz;

int main()
{
  // Create two 4x4 arrays.  We want them to look like matrices, so
  // we'll make the valid index range 1..4 (rather than 0..3 which is
  // the default).

  Range r(1,4);
  Array<float,2> A(r,r), B(r,r);

  // The first will be a Hilbert matrix:
  //
  // a   =   1
  //  ij   -----
  //       i+j-1
  //
  // Blitz++ provides a set of types { firstIndex, secondIndex, ... }
  // which act as placeholders for indices.  These can be used directly
  // in expressions.  For example, we can fill out the A matrix like this:

  firstIndex i;    // Placeholder for the first index
  secondIndex j;   // Placeholder for the second index

  A = 1.0 / (i+j-1);

  cout << "A = " << A << endl;
  // Now the A matrix has each element equal to a_ij = 1/(i+j-1).

  // The matrix B will be the permutation matrix
  //
  // [ 0 0 0 1 ]
  // [ 0 0 1 0 ]
  // [ 0 1 0 0 ]
  // [ 1 0 0 0 ]
  //
  // Here are two ways of filling out B:

  B = (i == (5-j));         // Using an equation -- a bit cryptic

  cout << "B = " << B << endl;

  B = 0, 0, 0, 1,           // Using an initializer list
    0, 0, 1, 0,           
    0, 1, 0, 0,
    1, 0, 0, 0;

  cout << "B = " << B << endl;

  // Now some examples of tensor-like notation.

  Array<float,3> C(r,r,r);  // A three-dimensional array: 1..4, 1..4, 1..4

  thirdIndex k;             // Placeholder for the third index

  // This expression will set
  //
  // c    = a   * b
  //  ijk    ik    kj

  C = A(i,k) * B(k,j);

  // In real tensor notation, the repeated k index would imply a
  // contraction (or summation) along k.  In Blitz++, you must explicitly
  // indicate contractions using the sum(expr, index) function:

  Array<float,2> D(r,r);

  D = sum(A(i,k) * B(k,j), k);

  // The above expression computes the matrix product of A and B.

  cout << "D = " << D << endl;

  // Now let's fill out a two-dimensional array with a radially symmetric
  // decaying sinusoid.

  int N = 64;                   // Size of array: N x N
  Array<float,2> F(N,N);
  float midpoint = (N-1)/2.;
  int cycles = 3;
  float omega = 2.0 * M_PI * cycles / double(N);
  float tau = - 10.0 / N;

  F = cos(omega * sqrt(pow2(i-midpoint) + pow2(j-midpoint)))
    * exp(tau * sqrt(pow2(i-midpoint) + pow2(j-midpoint)));

  return 0;
}

