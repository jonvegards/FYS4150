#include <iostream>
#include <vector>
#include <cmath>
#include <QElapsedTimer>
#include "planet.h"
#include "eulercromer.h"
using namespace std;

int main()
{
    vector<Planet*> planets;
    double sunMass = 1.0;
    double earthMass = 3.0024584e-6;
    Planet *earth = new Planet(vec3(1,0,0), vec3(0, 2*M_PI, 0), earthMass);
    Planet *sun = new Planet(vec3(0,0,0), vec3(0, 0, 0), sunMass);
    planets.push_back(earth);
    planets.push_back(sun);
    for(int i=0; i<100; i++) {
        double mass = rand()/RAND_MAX*1e-6;
        Planet *object = new Planet(vec3(), vec3(), mass);
        object->position.randomUniform(-2,2);
        object->velocity.randomUniform(-2,2);
        planets.push_back(object);
    }

    EulerCromer integrator;
    QElapsedTimer elapsedTimer;
    elapsedTimer.start();
    for(int t=0; t<1000; t++) {
        cout << t << endl;
        integrator.tick(planets, 1e-5);
    }
    cout << "Simulation finished after " << elapsedTimer.elapsed()/1000. << " seconds.";
    return 0;
}

// 0.738
