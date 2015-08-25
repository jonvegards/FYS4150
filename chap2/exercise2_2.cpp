/* Exercise 2.2 */
#include <stdlib.h> /* atof function atof = ascii to float*/
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function stdi = standard input output*/
#include <iostream>
using namespace std;
int main (int argc, char* argv[])
{
  	double N, i, j, k, sum1, sum2;        /* declare variables 64 bits*/
  	//N = atof(argv[1]);  /* convert the text argv[1] to double */
    double hurra[10], hurrb[10];
    int n, Q=5;
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
    double diff[10];
    for (n=0; n<=Q; n++){
      //cout << hurra[n] << endl;
      //cout << hurrb[n] << endl;
      // Calculating the logarithm of the difference
      diff[n] = log( (hurra[n]-hurrb[n])/hurrb[n] );
      cout << diff[n] << endl;
    }
    double final = sum2 - sum1;
    cout << final << endl;
  	return 0;           /* success execution of the program */
}