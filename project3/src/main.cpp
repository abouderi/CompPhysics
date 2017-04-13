#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"

using namespace std;

void PrintInitialValues(int, double, double, double *, double *, int);
void PrintFinalValues(int, double *, double *);

void PrintFinalValues(int Dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<Dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<Dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}

void PrintInitialValues(int Dimension,double TimeStep, double FinalTime,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << TimeStep << "; final time = " << FinalTime << "; integration points = " << N << endl;

    cout << "Initial position = ";
    for(int j=0;j<Dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<Dimension;j++) cout << v_initial[j] << " ";
    cout << endl;

}

int main()
{
    int IntegrationPoints;  // No. of integration points
    double FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

        cout << "Solar System" << endl;
        Dimension = 3;

        IntegrationPoints = 10000;
        FinalTime = 10000.; // units: years

        double TimeStep = FinalTime/((double) IntegrationPoints);
        double x[3],v[3];  // positions and velocities
        //double xJ[3],vJ[3]; // Jupiter positions and velocities (not necessary, but nice if you want to see them)
        // Initial Earth position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        // When adding planets, place them all in a line on the x-axis
        // Always place the Sun last. Idk why, it just works.
        // NOTE (PART C): Try sqrt(GM/r) for speed producing a circular orbit
        // GM = 4*M_PI*M_PI -> v = 2*M_PI/sqrt(r) = 2*M_PI
        // NOTE (PART D): Try sqrt(2GM/r) = sqrt(2)*2*M_PI for escape speed
        // NOTE (PART E): Simply change Jupiter's mass to see how it affects Earth's orbit
        // (mass,x,y,z,vx,vy,vz), Vel ~ 2*M_PI/sqrt(r)

        // contains distance (in AU) from sun for each body, in order of planet distance
        double r[10] = {0.0, 0.39, 0.72, 1.0, 1.52, 5.20, 9.54, 19.19, 30.06, 39.53};
        // contains masses of each body as a fraction of the Sun's mass, in order of distance
        double mass[10] = {1.0, 0.0000001659, 0.00000246, 0.000003, 0.0000003318, 0.0009546, 0.0002765, 0.00004424, 0.00005178, 0.000000006586};

        // NOTE (PART F): Center of Mass calculations for the whole system:
        double CoM = 0;
        double totalMass = 0;
        // Calculates center of mass
        for(int i = 0; i < 10; i++){
            CoM += r[i]*mass[i];
            totalMass += mass[i];
        }
        CoM = CoM/totalMass;
        // New positions relative to the center of mass of the solar system
        double r_new[10];
        for(int i = 0; i < 10; i++){
            r_new[i] = r[i] - CoM;
        }
        // Calculates new sun velocity to balance the momentum of the entire system
        // In the y-direction
        double sun_velocity = 0;
        for(int i = 1; i < 10; i++){
            sun_velocity -= 2*M_PI*(mass[i]/sqrt(r[i]));
        }
        // To do the new center of mass simulation (part F), simply replace r with r_new, and give the sun its sun_velocity
        planet Earth(mass[3],r[3],0,0.0,0,2*M_PI,0);
       // planet Jupiter(mass[5],r[5],0,0,0,2*M_PI/sqrt(r[5]),0);
        //planet Mercury(mass[1],r[1],0,0,0,2*M_PI/sqrt(r[1]),0);
        //planet Venus(mass[2],r[2],0,0,0,2*M_PI/sqrt(r[2]),0);
        //planet Mars(mass[4],r[4],0,0,0,2*M_PI/sqrt(r[4]),0);
        //planet Saturn(mass[6],r[6],0,0,0,2*M_PI/sqrt(r[6]),0);
        //planet Uranus(mass[7],r[7],0,0,0,2*M_PI/sqrt(r[7]),0);
        //planet Neptune(mass[8],r[8],0,0,0,2*M_PI/sqrt(r[8]),0);
        //planet Pluto(mass[9],r[9],0,0,0,2*M_PI/sqrt(r[9]),0);
        planet Sun(mass[0],r[0],0.,0.,0.,0.,0.);

        solver binary_vv(5.0);

        // Replicate these lines for each planet you add.
        binary_vv.add(Earth);
        //binary_vv.add(Jupiter);
        //binary_vv.add(Mercury);
        //binary_vv.add(Venus);
        //binary_vv.add(Mars);
        //binary_vv.add(Saturn);
        //binary_vv.add(Uranus);
        //binary_vv.add(Neptune);
        //binary_vv.add(Pluto);
        binary_vv.add(Sun);

        // Places initial Earth positions and velocities in their arrays
        for(int i = 0; i < Dimension; i++){
            x[i] = Earth.position[i];
            v[i] = Earth.velocity[i];
        }
 /*
        for(int i = 0; i < Dimension; i++){
            xJ[i] = Jupiter.position[i];
            vJ[i] = Jupiter.velocity[i];
        }
 */
        cout << "Initial Earth values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);
        cout << endl;
 /*
        cout << "Initial Jupiter values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,xJ,vJ,IntegrationPoints);
        cout << endl;
 */
        // Uncomment the second pair of lines to use the Euler algorithm
        cout << "Velocity Verlet results for the system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,1,0.);
        cout << "Velocity Euler results for the Sun-Earth system:" << endl;
        binary_vv.VelocityEuler(Dimension,IntegrationPoints,FinalTime,1,0.);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }
 /*
        for(int j = 0; j < Dimension;j++){
            xJ[j] = binary_vv.all_planets[1].position[j];
            vJ[j] = binary_vv.all_planets[1].velocity[j];
        }
 */
        PrintFinalValues(Dimension,x,v);
        cout << endl;
        //cout << "Jupiter positions and velocities: " << endl;
        //PrintFinalValues(Dimension,xJ,vJ);

    return 0;
}


// CENTER OF MASS SYSTEM
/*
int main()
{
    int IntegrationPoints;  // No. of integration points
    double FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

        cout << "Earth-Sun binary system" << endl;
        Dimension = 3;

        IntegrationPoints = 10000;
        FinalTime = 50.;
        // time_step = FinalTime/IntegrationPoints

        double TimeStep = FinalTime/((double) IntegrationPoints);
        double x[3],v[3];  // positions and velocities
        double xJ[3],vJ[3]; // Jupiter positions and velocities
        // initial position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        planet planet1(0.000003,0.9950378,0,0.0,0,2*M_PI,0); // Earth: (mass,x,y,z,vx,vy,vz)
        planet planet2(0.0009546,5.195,0,0,0,2*M_PI/sqrt(5.20),0); // Jupiter //vel ~ 2*M_PI/sqrt(5.20) ~ 2.7553590
        planet planet3(1.,-0.0049622,0.,0.,0.,-2*M_PI*(0.000003+0.0009546/sqrt(5.20)),0.); // Sun: (mass,x,y,z,vx,vy,vz)
        //NOTE (PART C): Try sqrt(GM/r) for speed producing a circular orbit
       // GM = 4*M_PI*M_PI -> v = 2*M_PI/sqrt(r) = 2*M_PI
       // NOTE (PART D): Try sqrt(2GM/r) = sqrt(2)*2*M_PI ~ 8.8857659 for escape speed
        solver binary_vv(5.0);
        binary_vv.add(planet1);
        binary_vv.add(planet2);
        binary_vv.add(planet3);

        //cout << "pos: " << x << endl;
        for(int i = 0; i < Dimension; i++){
            x[i] = planet1.position[i];
            v[i] = planet1.velocity[i];
        }

        for(int i = 0; i < Dimension; i++){
            xJ[i] = planet2.position[i];
            vJ[i] = planet2.velocity[i];
        }
        cout << "Initial Earth values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);
        cout << endl;

        cout << "Initial Jupiter values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,xJ,vJ,IntegrationPoints);
        cout << endl;

        cout << "Velocity Verlet results for the Sun-Jupiter system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,1,0.);
        //cout << "Velocity Euler results for the Sun-Earth system:" << endl;
        //binary_vv.VelocityEuler(Dimension,IntegrationPoints,FinalTime,1,0.);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }

        for(int j = 0; j < Dimension;j++){
            xJ[j] = binary_vv.all_planets[1].position[j];
            vJ[j] = binary_vv.all_planets[1].velocity[j];
        }

        PrintFinalValues(Dimension,x,v);
        cout << endl;
        cout << "Jupiter positions and velocities: " << endl;
        PrintFinalValues(Dimension,xJ,vJ);

    return 0;
}
*/
