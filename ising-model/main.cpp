/*
Program to solve the two-dimensional Ising model
The coupling constant J = 1
Boltzmann’s constant = 1, temperature has thus dimension energy
Metropolis sampling is used. Periodic boundary conditions.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <ctime>
using namespace std;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

// Function to read in data from screen
void read_input(int& n_spins, int& mcs, double& initial_temp, double& final_temp, double& temp_step){
    cout << "n_spins:";
    cin >> n_spins;

    cout << "mcs:";
    cin >> mcs;

    cout << "Initial Temp:";
    cin >> initial_temp;

    cout << "Final Temp:";
    cin >> final_temp;

    cout << "Temperature Step Size:";
    cin >> temp_step;
};


// Function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M)
{
	srand ( time(NULL) );





	
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){

	
	 int random = rand( )%2;
	double exponent = pow ( -1, rand( )%2);
	//cout << "Random: " << 1*exponent << endl;

	//cout << "This is the random number: " << endl;
	//cout << random << endl;
        spin_matrix[y][x] = exponent; // spin orientation for the ground state
        M += (double) spin_matrix[y][x];


        }
    }





    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        E -= (double) spin_matrix[y][x]*
        (spin_matrix[periodic(y,n_spins,-1)][x] +
        spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise

// The Metropolis algorithm
void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, int&matrixcounter, int&energycounter)
{

int babymatrixcounter = 0;


	
    // loop over all spins
        for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (ran1(&idum)*(double)n_spins);
            int iy = (int) (ran1(&idum)*(double)n_spins);

            int deltaE = 2*spin_matrix[iy][ix]*
            (spin_matrix[iy][periodic(ix,n_spins,-1)]+
            spin_matrix[periodic(iy,n_spins,-1)][ix] +
            spin_matrix[iy][periodic(ix,n_spins,1)] +
            spin_matrix[periodic(iy,n_spins,1)][ix]);


            // Here we perform the Metropolis test
                if ( ran1(&idum) <= w[deltaE+8] ) {
		babymatrixcounter = babymatrixcounter +1;		
                spin_matrix[iy][ix] *= -1;
                // flip one spin and accept new spin config
                // update energy and magnetization
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;


            }
	//cout << "Delta E: " << deltaE << endl;



        }
    }
	


	if(babymatrixcounter>0){
	matrixcounter = matrixcounter +1;
	//cout << matrixcounter << endl;
	}


} // end of Metropolis sampling over spins

void output(int n_spins, int mcs, double temp, double *average){
    cout << "n_spins:" << n_spins << endl;
    cout << "mcs:" << mcs << endl;
    cout << "temp:" << temp << endl;
    cout << "average:" << *average << endl;
};




// main program
int main(int argc, char* argv[])
{
	int matrixcounter;
	int energycounter;
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] <<
        "read also output file on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }
ofile.open(outfilename);

// Read in initial values such as size of lattice, temp and cycles
read_input(n_spins, mcs, initial_temp, final_temp, temp_step);

spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));

idum = -1; // random starting point

double integercounter = 1;
for (double temp = initial_temp; temp <= final_temp; temp+=temp_step){
    // initialise energy and magnetization
    matrixcounter = 0;
    energycounter = 0;
    E = M = 0.;

    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    initialize(n_spins, temp, spin_matrix, E, M);

double normal = mcs;
double spins = n_spins*n_spins;



//cout << "Number of spins: " << n_spins*n_spins << endl;








    // start Monte Carlo computation


    for (int cycles = 1; cycles <= mcs; cycles++){

	if(cycles == 1) int matrixcounter = 0;
	if(cycles == 1) int energycounter = 0;
	
        Metropolis(n_spins, idum, spin_matrix, E, M, w, matrixcounter, energycounter);
        // update expectation values
	

        average[0] += E/normal; average[1] += E*E/normal;
        average[2] += M/normal; average[3] += M*M/normal; average[4] += fabs(M/normal);
	//cout << "Energy: " << E << endl;
	

    } 




//define kb
double kb = 1;

double specificheat = 1/(kb*temp*temp)*(average[1]-average[0]*average[0]);
double chi = 1/(kb*temp)*(average[3]-average[4]*average[4]);


// print results
   // output(n_spins, mcs, temp, average);
/*
if(integercounter == 1) ofile << "Temperature " << "," << " Average Energy " << "," << " Average Abs Magnetization " << "," << " Specific Heat " << "," << " Susceptibility" << endl;
ofile << temp << "," << average[0] << "," << average[4] << "," << specificheat << "," << chi << endl;
*/
if(integercounter == 1) ofile << "#####################################" << endl;
ofile << "This is for run number " << integercounter << "!" << endl;
ofile << "Number of Spins: " << n_spins << endl;
ofile << "Monte Carlo Simulation Total: " << mcs << endl;
ofile << "Temperature: " << temp << endl;
ofile << "Average Energy: " << average[0] << endl;
ofile << "Average Energy Squared: " << average[1] << endl;
ofile << "Average Magnetization: " << average[2] << endl;
ofile << "Average Magnetization Squared: " << average[3] << endl;
ofile << "Absolute Average of Magnetization: " << average[4] << endl;
ofile << "" << endl;
ofile << "Specific Heat: " << specificheat << endl;
ofile << "Susceptibility: " << chi << endl;
ofile << "" << endl;
ofile << "There were " << matrixcounter << " different matrices produced." << endl;
ofile << "" << endl;
ofile << "The variance in E is: " << average[1]-average[0]*average[0] << endl;
ofile << "#####################################" << endl;


integercounter = integercounter +1;

}

cout << "Computation Finished" << endl;

free_matrix((void **) spin_matrix); // free memory
ofile.close(); // close output file
return 0;
}
