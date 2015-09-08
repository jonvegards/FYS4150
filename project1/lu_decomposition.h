#ifndef LU_DECOMPOSITION
#define LU_DECOMPOSITION

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>    // For printing max value in an array
#include <armadillo>    // How do I compile with this included?

double lu_decomposition(double *a, double *b, double *c, double *x, int n);

#endif // LU_DECOMPOSITION

