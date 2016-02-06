#include "planet.h"

Planet::Planet(vec3 initialPosition, vec3 initialVelocity, double massIn) :
    position(initialPosition),
    velocity(initialVelocity),
    mass(massIn)
{

}

