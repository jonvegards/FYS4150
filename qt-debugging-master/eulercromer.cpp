#include "eulercromer.h"
#include <cmath>

EulerCromer::EulerCromer()
{

}

void EulerCromer::calculateForces(vector<Planet*> planets) {
    for(Planet *planet : planets) {
        planet->force.setToZero();
    }

    double G = 4*M_PI*M_PI;
    for(int i=0; i<planets.size(); i++) {
        Planet *planet1 = planets[i];
        for(int j=i+1; j<planets.size(); j++) {
            Planet *planet2 = planets[j];
            double massProduct = planet1->mass*planet2->mass;
            vec3 &p1 = planet1->position;
            vec3 &p2 = planet2->position;

            double dx = p1.x() - p2.x();
            double dy = p1.y() - p2.y();
            double dz = p1.z() - p2.z();
            double dr2 = dx*dx + dy*dy + dz*dz;

            double force = G*massProduct/(dr2*sqrt(dr2));
            planet1->force[0] -= force*dx;
            planet1->force[1] -= force*dy;
            planet1->force[2] -= force*dz;

            planet2->force[0] += force*dx;
            planet2->force[1] += force*dy;
            planet2->force[2] += force*dz;


//            vec3 deltaR = (planet1->position-planet2->position);
//            double distance = deltaR.length();

//            double force = G*massProduct/pow(distance,3);

//            vec3 forceVector = force*deltaR;
//            planet1->force -= forceVector;
//            planet2->force += forceVector;
        }
    }
}

void EulerCromer::tick(vector<Planet*> planets, double dt) {
    calculateForces(planets);
    for(Planet *planet : planets) {
        planet->velocity += planet->force*dt;
        planet->position += planet->velocity*dt;
    }
}
