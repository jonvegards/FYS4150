/* Exercise 2.2 */
#include <stdlib.h> /* atof function atof = ascii to float*/
#include <math.h>   /* sine function */
#include <cmath>    /* absolute value */
#include <stdio.h>  /* printf function stdi = standard input output*/
#include <fstream>  /* to read/make files */
#include <iostream>
using namespace std;
int main (int argc, char* argv[])
{
  	double N, i, j, k, sum1, sum2;        /* declare variables 64 bits*/
  	//N = atof(argv[1]);  /* convert the text argv[1] to double */
    int n, Q=10;
    // Defining all necessary arrays
    double hurra[Q], hurrb[Q], diff[Q];
    int number[Q];

  	for (n=0; n<=Q; n++){
      N = pow(10,n);
      k = N;
      for (i=1; i<=N; i++){
        sum1 += 1/i;    // Sum up
        sum2 += 1/k;    // Sum down
        k -= 1;
      }
      hurra[n] = sum1;
      hurrb[n] = sum2;
    }

    // Making a file to send our data to
    ofstream file("exercise2_2.dat");
    file << "#x y" << endl;

    for (n=0; n<=Q; n++){
      // Calculating the logarithm of the difference
      number[n] = n;
      diff[n] = log( abs( (hurra[n]-hurrb[n])/hurrb[n] ) );
      cout << diff[n] << endl;
      file << number[n] << ' ' << diff[n] << endl;
    }
  	return 0;           /* success execution of the program */
}