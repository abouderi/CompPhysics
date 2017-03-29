#ifndef PLANET
#define PLANET
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include<iostream>


using std::vector;

class planet
{
public:
    double mass;
    double position[2];
    double velocity[2];
    double potential;
    double kinetic;
    
    
    
    planet();
    planet(double M, double x, double y, double vx, double vy);
    
    
    
    double distance(planet otherPlanet);
    double GravitationalForce(planet otherPlanet, double Gconst);
    double Acceleration(planet otherPlanet, double Gconst);
    double KineticEnergy();
    double PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon);
    
    
};

#endif