// Testfil i .cpp-stil for å sjå om giten er i form

// A comment line begins like this in C++ programs
// Standard ANSI-C++ include files
// using namespace std;
#include <iostream>  // input and output
#include <math.h>   // sine function
int main (int argc, char* argv[])
{
  // convert the text argv[1] to double using atof:
  double r = atof(argv[1]);
  double s = sin(r);
  std::cout << "Hei, verd! Har du det bra? sin(" << r << ")=" << s << '\n';
  // success
  return 0;
}