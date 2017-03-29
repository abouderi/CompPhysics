#include "planet.h"


planet::planet()
{
    mass = 1.;
    position[0] = 1.;
    position[1] = 0.0;
    //position[2] = 0.0;
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    //velocity[2] = 0.0;
    U = 0.0;
    KE = 0.0;
    
}

planet::planet(double M, double x, double y, double vx, double vy)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    velocity[0] = vx;
    velocity[1] = vy;
    U = 0.0;
    KE = 0.0;
    
}

double planet::distance(planet planetother)
{
    double x1,y1,x2,y2,xx,yy;
    
    x1 = this->position[0];
    y1 = this->position[1];
    
    x2 = planetother.position[0];
    y2 = planetother.position[1];
    
    xx = x1-x2;
    yy = y1-y2;
    
    return sqrt(xx*xx+yy*yy);
}

double planet::GravitationalForce(planet planetother, double G)
{ 
    double r = this->distance(planetother);
    
    if(r!=0) return G*this->mass*planetother.mass/(r*r);
    
    else return 0;
}

double planet::Acceleration(planet planetother,double G)
{
    double r = this->distance(planetother);
    if(r!=0) return this->GravitationalForce(planetother,G)/(this->mass*r);
    else return 0;
}

double planet::KE()
{
    double velocity2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]);
    return = 0.5*this->mass*velocity2;
}

double planet::U()
{
    if(epsilon==0) return -G*this->mass*planetother.mass/this->distance(planetother);
    else return(G*this->mass*planetother.mass*epsilon)*(atan(this->distance(planetother)/epsilon)-0.5*M_PI);
}
