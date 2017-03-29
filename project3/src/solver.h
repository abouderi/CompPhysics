#ifndef SOLVER
#define SOLVER
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>



using std::vector;

class solver
{
public:
    friend class planet;
    double radius,TotalMass,G;
    int TotalPlanets;
    vector<planet> AllPlanets;
    double KETotal;
    double UTotal;
    
    
    
    solver();
    solver(double radi);
    
    void add(planet newplanet);
    void addM(planet newplanet);
    void GraviationalConstant();
    void PrintPosition(std::ofstream &output, int dimension, double time, int number);
    void print_energy(std::ofstream &output, double time, double epsilon);
    void VelocityVerlet(int dimension, int integration_points, double final_time, int print_number, double epsilon);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double epsilon);
    void GravitationalForce_RK(double x_rel, double y_rel, double &Fx, double &Fy, double mass1, double mass2);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    double EnergyLoss();
    bool Bound(planet OnePlanet);


};


#endif