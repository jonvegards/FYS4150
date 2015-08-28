// Exercise 3.1
// Program that computes the first derivative to
// f(x) = arctan(x) for x = sqrt 2, steplength h
// Exact answer is 1/3.

/*
Using the two formulae
f_2c '(x) = (f(x+h) - f(x)) / (h),
f_3c '(x) = (f(x+h) - f(x-h)) / (2h)

x = sqrt 2
h = 0.1, this we'll improve afterwards
*/

#include <cmath> 		// Math functions
#include <iostream>		// Dunno
using namespace std;

// Allocating variables

double x;
double step_length[] = { 1, 0.1, 0.01, 0.001, 0.0001, 0.00001};
double computed_derivative1[sizeof(step_length)];
double computed_derivative2[sizeof(step_length)];

// Defining function for three-point formulae
// which takes two variables
void first_derivative1(double, double *, double *);
void first_derivative2(double, double *, double *);
void print_result(double *, double *, double);

int main()
{
    // Defining x-value
    x = sqrt(2);

    // Calculating the first derivative

    first_derivative1(x, computed_derivative1, step_length);
    first_derivative2(x, computed_derivative2, step_length);
    print_result(computed_derivative1, step_length, x);
    //print_result(computed_derivative2, step_length, x);
    return 0;
}

// Defining function to calculate the derivative of arctan(x)
void first_derivative1(double q, double *computed_derivative1, double *step_length){
    int counter;
    double w;

    for (counter=0; counter < sizeof(step_length); counter++){
        w = step_length[counter];
        computed_derivative1[counter] = ( atan(q+w) - atan(q) ) / (w);
    }
    return;
}

// Defining function to calculate the derivative of arctan(x) with higher precision
void first_derivative2(double e, double *computed_derivative2, double *step_length){
    int counter;
    double r;

    for (counter=0; counter < sizeof(step_length); counter++){
        r = step_length[counter];
        computed_derivative2[counter] = ( atan(e+r) - atan(e-r) ) / (2*r);
    return;
    }
}

void print_result(double *computed_derivative, double *step_length, double q){
    // Function for printing result

    FILE *output_file;
    output_file = fopen("out.txt", "w") ;

    int i;
    for (i=0; i<sizeof(computed_derivative); i++){
        fprintf(output_file, "%12.5E %12.5E \n",
                 log10(step_length[i]),log10(fabs(computed_derivative[i]-atan(q))/atan(q)));
    }
    fclose (output_file);
}
