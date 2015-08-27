/* Exercise 2.2 */
#include <stdlib.h> /* atof function atof = ascii to float*/
#include <math.h>   /* sine function */
#include <cmath>    /* absolute value */
#include <stdio.h>  /* printf function stdi = standard input output*/
#include <matplotpp.h> /* plotting */
#include <iostream>
using namespace std;
int main (int argc, char* argv[])
{
  	double N, i, j, k, sum1, sum2;        /* declare variables 64 bits*/
  	//N = atof(argv[1]);  /* convert the text argv[1] to double */
    int n, Q=5;
    double hurra[Q], hurrb[Q];
    int number [] = { 0, 1, 2, 3, 4 }; //, 5, 6, 7, 8, 9 };
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
  	//printf("Sum up: %2.8f" ,sum1); /* Print result */
    //cout << endl; 
  	//printf("Sum up: %2.8f" ,sum2); /* Print result */
  	//cout << endl;
    double diff[Q];
    for (n=0; n<=Q; n++){
      //cout << hurra[n] << endl;
      //cout << hurrb[n] << endl;
      // Calculating the logarithm of the difference
      diff[n] = log( abs( (hurra[n]-hurrb[n])/hurrb[n] ) );
      cout << diff[n] << endl;
    }
    double final = sum2 - sum1;
    cout << final << endl;

    // Plotting the results
    plotxy()

  	return 0;           /* success execution of the program */
}