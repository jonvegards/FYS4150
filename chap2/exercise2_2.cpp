/* Exercise 2.2 */
#include <stdlib.h> /* atof function atof = ascii to float*/
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function stdi = standard input output*/
#include <iostream>
using namespace std;
int main (int argc, char* argv[])
{
  	double N, n, i, j, sum1, sum2;        /* declare variables 64 bits*/
  	//N = atof(argv[1]);  /* convert the text argv[1] to double */
  	N = pow(10,10);
	n = N;
  	// Sum up
  	for (i=1; i<=N; i++)
  	{
  		sum1 += 1 / i;
  		sum2 += 1/n; // Sum down
  		n -= 1;
  	}
  	printf("Sum up: %2.8f" ,sum1); /* Print result */
    cout << endl; 
  	printf("Sum up: %2.8f" ,sum2); /* Print result */
  	cout << endl;
    double final = sum2 - sum1;
    cout << final;
    cout << endl;
  	return 0;           /* success execution of the program */
}