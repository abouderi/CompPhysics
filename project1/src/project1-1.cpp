#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iomanip> 
#include <string>
#include "time.h"   //  you have to include the time.h header


using namespace std;  // Why namespace? 

// object for output files, for input files use ifstream
ofstream ofile;
// Functions used, why is this an inline function?
inline double ff(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program
int main(int argc, char *argv[]){
  int exponent; 
  string filename; //  This is a useful way to define a text string, alternatively you can use characters
    // We read also the basic name for the output file and the highest power of 10^n we want
  
  clock_t start, finish;  //  declare start and final time
    start = clock();
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
    }

      //Many lines of comments start like this
    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);

      //   To be filled in during lab session
      double h = 1.0/((double) n); //make sure we don't have integer division
      double hh = h*h;
      
      //need to have flexible (dynamic memory allocation)
      double *d = new double[n+1];  //diagonal matrix elements
      double *f = new double[n+1];  //Right Hand Side
      //double *e = new double[n+1];  //non-diagonal matrix elements
      double *x = new double[n+1];  //x-value
      double *u = new double[n+1];  //solution
      
      u[0] = 0.0;  u[n] = 0.0;  //u[0] = u[n] = 0.0;
      
      d[0] = d[n] = 2.0;
      
      //
      for (int i=0; i <= n; i++){
          x[i] = i*h;
          f[i] = hh*ff(i*h);
          
      }
      
      
      //set up d[i] where d[i]=(1+i)/i
      for (int i=1; i<n; i++){
          d[i] = (i+1.0)/((double) i);
      }
      
      //Forward Substitution
      for (int i=2; i<n; i++){
          f[i] = f[i] + f[i-1]/d[i-1];
      }
      //Same thing as above
      //for (int i=2; i<n; i++) f[i] += f[i-1]/d[i-1];
   
      
      
      //Backward Substitution
      
      u[n-1] = f[n-1]/d[n-1];
      
      for (int i = n-2; i>0; i--) u[i] = (f[i]+u[i+1])/d[i];
      
      



      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      
      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-u[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << u[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] d; delete [] f; delete [] u;
     
  
      
      finish = clock();
      cout << "Time:" << endl;
    cout << ( (finish - start)/CLOCKS_PER_SEC ) << endl;
      
      
      
    }
      
    return 0;   // returns zero upon success
}
