#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "time.h"   //  you have to include the time.h header


using namespace std;

#define ZERO    1.0E-15

/* function declarations */
double ** AllocateMatrix(int, int);
void DeallocateMatrix(double **, int, int); 
void MatrixInverse(double **, int);
void WriteMatrix(double **, int);
void MatrixMultiplication(double **, double **, int);
void LUDecomposition(double **, int, int *);
void LUBackwardSubstitution(double **, int, int *, double *);



int main(int argc, char *argv[]){
    int n = atoi(argv[1]);
    
    double **A = AllocateMatrix(n, n);
    double **B = AllocateMatrix(n ,n);
    
    double invers_period = 1./RAND_MAX;
    
    
    clock_t start, finish;  //  declare start and final time
    start = clock();
    
    srand(time(NULL));
    // Setting general square matrices with random matrix elements
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            double x= 0.0;
            A[i][j] = x;
            
          
            B[i][j] = A[i][j];
            
           
        }
    }
    
    
    

cout << "Initial A:" << endl;
WriteMatrix(A, n);

MatrixInverse(B, n);

cout << "Inverse A:" << endl;
WriteMatrix(B, n);
 
cout << "We should have an identity matrix:" << endl;
MatrixMultiplication(A,B,n);

finish = clock();
cout << "Time:" << endl;
 cout <<   ( (finish - start)/CLOCKS_PER_SEC ) << endl;

return 0;

}



void MatrixInverse(double **A, int n)
{        
  // allocate space in memory
  int *indx;  
  double *column;
  indx = new int[n];
  column  = new double[n];
  double **Y    = AllocateMatrix(n,n);
  // Perform the LU decomposition 
  LUDecomposition(A, n, indx);   // LU decompose  a[][] 
  cout << "LU decomposed matrix  A:" << endl;
  WriteMatrix(A,n);
  // find inverse of a[][] by columns 
  for(int j = 0; j < n; j++) {
    // initialize right-side of linear equations 
    for(int i = 0; i < n; i++) column[i] = 0.0;
    column[j] = 1.0;
    LUBackwardSubstitution(A, n, indx, column);
    // save result in y[][] 
    for(int i = 0; i < n; i++) Y[i][j] = column[i];
  } 
  // return the inverse matrix in A[][] 
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) A[i][j] = Y[i][j];
  } 
  DeallocateMatrix(Y, n, n);     // release local memory 
  delete [] column;
  delete []indx;

} 











double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}

// Free memory

void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}

// Write out a given matrix
void WriteMatrix(double ** Matrix, int n){
  for(int i=0;i < n;i++){
    cout << endl;
     for (int j=0 ; j < n;j++){
        printf("  A[%2d][%2d] = %12.4E",i, j, Matrix[i][j]);
     }
  }
    cout << endl;
}

// Straightforward matrix-matrix multiplication

void MatrixMultiplication(double ** a, double **b, int n){
  double **c = AllocateMatrix(n, n);
  for(int i=0;i < n;i++){
     for (int j=0 ; j < n;j++){
       double sum = 0.0;
       for (int k = 0; k < n; k++) sum += a[i][k]*b[k][j];
       c[i][j] = sum;
     }
  }
  WriteMatrix(c,n);
}

/*
    The function
    void LUDecomposition(double **a, int n, int *indx)
    takes as input a two-dimensional matrix a[][] of dimension n and
    replaces it by the LU decomposition of a rowwise permutation of
    itself. The results is stored in a[][]
    The vector
    indx[] records the row permutation effected by the partial pivoting.
*/

void LUDecomposition(double **a, int n, int *indx)
{
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new double [n];
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a[i][j])) > big) big = temp;
      }
      if(big == ZERO) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }               
      vv[i] = 1.0/big;        
   } 

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a[i][j];    
	 for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
	 a[i][j] = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a[i][j];
	 for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
	 a[i][j] = sum;
	 if((dum = vv[i]*fabs(sum)) >= big) {
  	    big = dum;
	    imax = i;
	 }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
	    dum        = a[imax][k];
	    a[imax][k] = a[j][k];
	    a[j][k]    = dum;
	 }
	 vv[imax] = vv[j];         // also interchange scaling factor 
      }
      indx[j] = imax;
      if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;
      if(j < (n - 1)) {                   // divide by pivot element 
         dum = 1.0/a[j][j];
	 for(i=j+1;i < n; i++) a[i][j] *= dum;
      }
   } // end j-loop over columns
   delete [] vv;   // release local memory
}
 

/*
     The function 
       void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
     solves the set of linear equations A X = B of dimension n.
     a[][] is input, not as the matrix A[][] but rather as 
     its LU decomposition, indx[] is input as the permutation vector returned by 
     ludcmp(). b[] is input as the right-hand side vector B,
     The solution X is returned in B. The input data a[][],
     n and indx[] are not modified. This routine takes into 
     account the possibility that b[] will begin with many
     zero elements, so it is efficient for use in matrix
     inversion.
*/

void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
{
   int        i, ii = -1, ip, j;
   double     sum;

   for(i = 0; i< n; i++) {
      ip    = indx[i];
      sum   = b[ip];
      b[ip] = b[i];
      if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
      else if(sum) ii = i;
      b[i] = sum;
   }
   for(i = n - 1; i >= 0; i--) {
      sum = b[i];
      for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
      b[i] = sum/a[i][i];
   }
}