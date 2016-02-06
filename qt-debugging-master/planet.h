#ifndef PLANET_H
#define PLANET_H
#include "vec3.h"

class Planet
{
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass = 0;
    Planet(vec3 initialPosition, vec3 initialVelocity, double mass);
};

#endif // PLANET_H
