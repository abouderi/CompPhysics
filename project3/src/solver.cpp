#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"

/*
  double radius,TotalMass,G;
    int TotalPlanets;
    vector<planet> AllPlanets;
    double KETotal;
    double UTotal;
 */

solver::solver()
{
    TotalPlanets = 0;
    radius = 100;
    TotalMass = 0;
    G = 4*M_PI*M_PI;
    KETotal = 0;
    UTotal = 0;
    
}

solver::solver(double radi)
{
   TotalPlanets = 0;
   radius = radi;
   TotalMass = 0;
   G = 4*M_PI*M_PI;
   KETotal = 0;
   UTotal = 0;
   
}

void solver::add(planet newplanet)
{
    TotalPlanets +=1;
    AllPlanets.pushback(newplanet);
    
}

void solver::GravitationalConstant()
{
     G = (4*M_PI*M_PI/32)*radius*radius*radius/TotalMass;
}

void solver::PrintPosition(std::ofstream &output, int dimension, double time, int number)
{
    if(dimension>2 || dimension <= 0) dimension = 2;
    else{
        for(int i=0;i<number;i++){
            planet &current = AllPlanets[i];
            output << time << "\t" << i+1 << "\t" << current.mass;
            for(int j=0;j<dimension;j++) output << "\t" << current.position[j];
            for(int j=0;j<dimension;j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
    }
}


void solver::print_energy(std::ofstream &output, double time,double epsilon)
{   // Writes energies to a file "output"

    this->KineticEnergySystem();
    this->PotentialEnergySystem(epsilon);
    for(int nr=0;nr<TotalPlanets;nr++){
        planet &Current = AllPlanets[nr];
        output << time << "\t" << nr << "\t";
        output << Current.KE << "\t" << Current.U << std::endl;
    }
}

void solver::VelocityVerlet(int dimension, int integration_points, double final_time, int print_number, double epsilon)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/


    double TimeStep = FinalTime/((double) IntegrationPoints);
    double time = 0.0;
    double loss = 0.; //Energy Loss
    int lostPlanets[IntegrationPoints];


    char *filename = new char[1000];
    char *filenameE = new char[1000];
    char *filenameB = new char[1000];
    char *filenameLost = new char[1000];
        sprintf(filename, "PlanetsVV_%d_%.3f.txt",total_planets,time_step); 
        sprintf(filenameE, "PlanetsVV_energy_%d_%.3f.txt",total_planets,time_step);
        sprintf(filenameB,"Planetsbound_%d_%.3f.txt",total_planets,time_step);
        sprintf(filenameLost,"Planetslost_%d_%.3f.txt",total_planets,time_step);
    std::ofstream output_file(filename);
    std::ofstream output_energy(filenameE);
    std::ofstream output_bound(filenameB);
    std::ofstream output_lost(filenameLost);


    double **Acceleration = setup_matrix(TotalPlanets,3);
    double **Acceleration_new = setup_matrix(TotalPlanets,3);


    double Fx,Fy,Fxnew,Fynew;


    PrintPosition(output_file,dimension,time,print_number);
    print_energy(output_energy,time,epsilon);

    int n = 0;
    lostPlanets[n] = 0;
    output_lost << time << "\t" << lostPlanets[n] << std::endl;
    n+=1;

    clock_t planet_VV,finish_VV;
    planet_VV = clock();


    time += time_step;
    while(time < final_time){
        lostPlanets[n] = 0;


        for(int nr1=0; nr1<TotalPlanets; nr1++){
            planet &current = AllPlanets[nr1];

            Fx = Fy = Fxnew = Fynew = 0.0;


                for(int nr2=nr1+1; nr2<TotalPlanets; nr2++){
                    planet &other = AllPlanets[nr2];
                    GravitationalForce(current,other,Fx,Fy,epsilon);
                }

            Acceleration[nr1][0] = Fx/current.mass;
            Acceleration[nr1][1] = Fy/current.mass;
            //Acceleration[nr1][2] = Fz/current.mass;


            for(int j=0; j<dimension; j++) {
                current.position[j] += current.velocity[j]*TimeStep + 0.5*TimeStep*TimeStep*Acceleration[nr1][j];
            }


                for(int nr2=nr1+1; nr2<TotalPlanets; nr2++){
                    planet &other = AllPlanets[nr2];
                    GravitationalForce(current,other,Fxnew,Fynew,epsilon);
                }


            Acceleration_new[nr1][0] = Fxnew/current.mass;
            Acceleration_new[nr1][1] = Fynew/current.mass;
           // acceleration_new[nr1][2] = Fznew/current.mass;


            for(int j=0; j<dimension; j++) current.velocity[j] += 0.5*TimeStep*(Acceleration[nr1][j] + Acceleration_new[nr1][j]);
        }



        pPrintPosition(output_file,dimension,time,print_number);
        print_energy(output_energy,time,epsilon);

        loss += EnergyLoss();

        for(int nr=0;nr<TotalPlanets;nr++){
            planet &Current = AllPlanets[nr];
            if(!(this->Bound(Current))){
                lostPlanets[n] += 1;
            }
        }
        output_lost << time << "\t" << lostPlanets[n] << std::endl;
        n += 1;
        time += time_step;
    }
  
    finish_VV = clock();
    std::cout << "Total time = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl;
    std::cout << "One time step = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC)/IntegrationPoints << " seconds" << std::endl; 

    //loss = EnergyLoss();
    std::cout << "Total energyloss due to unbound planets: " << loss << std::endl;

    double boundPlanets = 0;
    for(int nr=0;nr<TotalPlanets;nr++){
        planet &Current = AllPlanets[nr];
        if(this->Bound(Current)){
            output_bound << nr << std::endl;
            boundPlanets += 1;
        }
    }
    std::cout << "There are " << boundPlanets << " bound planets at the end of the run" << std::endl;

    // Close files
    output_file.close();
    output_energy.close();
    output_bound.close();
    output_lost.close();

    // Clear memory
    delete_matrix(acceleration);
    delete_matrix(acceleration_new);
}

double ** solver::setup_matrix(int height, int width)
{
    double **matrix;
    matrix = new double*[height];
    
    
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];
        
    for(int i=0;i<height;i++)
    {
       for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        } 
    }
        
       return matrix; 
        
}

void solver::delete_matrix(double **matrix)
{
    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void solver::GravitationalForce(planet &current,planet &other,double &Fx,double &Fy,double epsilon)
{   

    double relative_distance[2];

    for(int j = 0; j < 2; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    double smoothing = epsilon*epsilon*epsilon;


    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r) + smoothing);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r) + smoothing);
    //Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r) + smoothing);
}

void solver::GravitationalForce_RK(double x_rel,double y_rel,double &Fx,double &Fy,double mass1,double mass2)
{   


    double r = sqrt(x_rel*x_rel + y_rel*y_rel);

   
    Fx -= this->G*mass1*mass2*x_rel/(r*r*r);
    Fy -= this->G*mass1*mass2*y_rel/(r*r*r);
    //Fz -= this->G*mass1*mass2*z_rel/(r*r*r);
}

void solver::KineticEnergySystem()
{
    KETotal = 0;
    for(int nr=0;nr<TotalPlanets;nr++){
        planet &Current = AllPlanets[nr];
        Current.kinetic = Current.KE();
    }
}

void solver::PotentialEnergySystem(double epsilon)
{
    UTotal = 0;
    for(int nr=0;nr<TotalPlanets;nr++){
        planet &Current = AllPlanets[nr];
        Current.potential = 0;
    }
    for(int nr1=0;nr1<TotalPlanets;nr1++){
        planet &Current = AllPlanets[nr1];
        for(int nr2=nr1+1;nr2<TotalPlanets;nr2++){
            planet &Other = AllPlanets[nr2];
            Current.potential += Current.U(Other,G,epsilon);
            Other.potential += Other.U(Current,G,epsilon);
        }
    }
}

bool solver::Bound(planet OnePlanet)
{
    return ((OnePlanet.kinetic + OnePlanet.potential) < 0.0);
}


double solver::EnergyLoss()
{
    bool bound;
    vector<int> indices;
    double EnergyLoss = 0;
    for(int nr=0;nr<TotalPlanets;nr++){
        planet &Current = AllPlanets[nr];
        bound = this->Bound(Current);
        if(!bound){
            indices.push_back(nr);
        }
    }
    for(int i=0;i<indices.size();i++) EnergyLoss += AllPlanets[indices[i]].KE();
    return EnergyLoss;
}