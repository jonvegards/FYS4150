#include <iostream>
#include "vec3.h"

using namespace std;

int main()
{
    cout << "Hello World!" << endl;
    vec3 myVector(1,2,3);
    myVector.print();
    cout << myVector.length() << endl;
    cout << myVector[0] << endl;
    cout << myVector.x() << endl;
    vec3 a(1,2,3);
    vec3 b(2,3,4);
    vec3 c = a + b;
    c.print();
    cout << (a+b).length() << endl;
    vec3 d = a + 2;
    d.print();

    return 0;
}
