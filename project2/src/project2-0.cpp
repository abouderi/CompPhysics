#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iomanip> 
#include <string>
#include "time.h" 
#include <armadillo>



using namespace std;
using namespace arma;


int main(int argc, char *argv[]){
    int n = atoi(argv[1]);
    
    vec V(n);
    
    
    
    V(n) = orth(V(n));
    
    
    
    
    
    
    
    
    return 1;
}