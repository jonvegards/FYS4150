// Trying to learn some function stuff bby
#include <iostream>
#include <stdlib.h> /* atof function */
using namespace std;

// First, we define an inline function
inline double MIN(double a,double b) {return (((a)<(b)) ? (a):(b));}

int main(int argc, char *argv[])
{
	/* Just calling on the inline function */
	cout << MIN(atof(argv[1]),atof(argv[2])) << endl;
	return 0;
}