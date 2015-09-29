// Testing code for project 2
#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

int main()
{
	vec r1 = zeros<vec>(5);
	vec r2 = randu<vec>(5);
	vec r = r1 - r2;
	r1.print();
	r2.print();
	r.print();
	return 0;
}