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
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

//-larmadillo -llapack -lblas

using namespace std;
using namespace arma;

ofstream ofile;

#define ZERO    1.0E-15



double offdiag(mat &A , int &p, int &q, int n){
     //int p, q;
         double max = 0.0;
         for(int i=0; i<n; ++i){
                for(int j=i+1; j<n; ++j)
                {
                    double aij = fabs(A(i,j));
                    if(aij>max){
                       max = aij; p=i; q=j;
                         }
                 }
             }
         return max;
}


void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n )
{
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
    
   
  }
  return;
} // end of function jacobi_rotate
   






int main(int argc, char *argv[]){
    int n = atoi(argv[1]);
   
    double rho_max = 7.0;
    double rho_min = 0.0;
    
    double h = (rho_max - rho_min)/n;
    
    //cout << h << endl;
    
    double hh = h*h;
    
    //cout << hh << endl;
    
    double hh_inverse = 1.0/hh;
    
    //cout << hh_inverse << endl;
    
    
    
    
    
    clock_t start, finish;  //  declare start and final time
    start = clock();
    
    srand(time(NULL));
    
    
    
    vec rho(n); vec V(n);
    
    double omegr = 5;
    
    //Defining Rho
    for(int i=0; i<n; i++){
        
        rho(i) = rho_min+(i+1)*h;
        V(i) = omegr*omegr*rho(i)*rho(i)+1/rho(i);
        
    }
    
    //cout << rho << endl;

    
    //cout << V << endl;
    
    mat A = zeros<mat>(n,n);
    
    // Setting the diagonal and off diagonal elements.
    for(int i=0; i<n-1; i++){
            A(i,i) = (2.0*hh_inverse)+V[i];
            A(i+1,i+1) = (2.0*hh_inverse)+V[i+1];
            A(i,i+1) = -1.0*hh_inverse;
            A(i+1,i) = -1.0*hh_inverse;
            
            
    }
    
        vec vectorprior(n);
        vec vectorprior2(n);
        
       for(int y=0; y<n-1; y++){
            //vectorprior(y) = A(y,0);
            if (A(y,0)<0.0){
                vectorprior(y) = 0.0;
            }else{
                vectorprior(y) = A(y,0);
            }
            if (A(y,1)<0.0){
                vectorprior2(y) = 0.0;
            }else{
                vectorprior2(y) = A(y,1);
            }
                
            //vectorprior2(y) = A(y,1);
        }
        
        cout << "Eigenvector One (Before Jacobian Trans): " << endl;
        cout << vectorprior << endl;
        cout << "Eigenvector Two (Before Jacobian Trans): " << endl;
        cout << vectorprior2 << endl;
        
        double c = dot(vectorprior,vectorprior2);
        
        double w = dot(vectorprior,vectorprior);
        
        cout << "Dot Product of Two Eigenvectors: " << c << endl;

        
        cout << "Dot Product of One Eigenvector: " << w << endl;
    
    cout << "Matrix A initial:" << endl;
    cout << A << endl;

    mat A1 = A;
    
    vec eigval;

    eig_sym(eigval, A1);
    
    cout << "Eigenvalues via Armadillo:" << endl;
    cout << eigval << endl;

    mat R = eye<mat>(n,n);
    
    
    double maxnondiag = 1;
    double maxiter = 100;
    double tolerance = 1.0E-10;
    int iterations = 0;
    while(maxnondiag > tolerance && iterations <= maxiter )
    {
 
        int p, q;
        maxnondiag = offdiag(A,p,q,n);
        //cout << maxnondiag << endl;
        //cout << p << q << endl;
        Jacobi_rotate(A,R,p,q,n);
        iterations++;
        
    
        

    
    }
 
        cout << "Matrix A after Jacobi Rotation:" << endl;
        cout << A << endl;
        
        vec vector(n);
        vec vector2(n);
        
        for(int y=0; y<n-1; y++){
            //vector(y) = A(y,0);
            //vector2(y) = A(y,1);
            
            if (A(y,0)<0.0){
                vector(y) = 0.0;
            }else{
                vector(y) = A(y,0);
            }
            if (A(y,1)<0.0){
                vector2(y) = 0.0;
            }else{
                vector2(y) = A(y,1);
            }
        }
        
        cout << "Eigenvector One: " << endl;
        cout << vector << endl;
        cout << "Eigenvector Two: " << endl;
        cout << vector2 << endl;
        
        double x = dot(vector,vector2);
        
        double q = dot(vector,vector);
        
        cout << "Dot Product of Two Eigenvectors (After Jacobian Trans): " << x << endl;

        
        cout << "Dot Product of One Eigenvector (After Jacobian Trans): " << q << endl;
        
        

        
 
    return 0;
    
}