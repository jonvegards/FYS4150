#ifndef EULERCROMER_H
#define EULERCROMER_H
#include <vector>
#include "planet.h"
using std::vector;

class EulerCromer
{
public:
    EulerCromer();
    void calculateForces(vector<Planet*> planets);
    void tick(vector<Planet*> planets, double dt);
};

#endif // EULERCROMER_H
