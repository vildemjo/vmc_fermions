// # include "matrix_inverse.h"


// std::vector<std::vector<double>> ludcmp(double **a, int n, int *indx, double *d)
// {
//    int      i, imax, j, k;
//    double   big, dum, sum, temp, *vv;

//   vv = new(nothrow) double [n];
//   if(!vv) {
//     printf("\n\nError in function ludcm():");
//     printf("\nNot enough memory for vv[%d]\n",n);
//     exit(1);
//   }

//    *d = 1.0;                              // no row interchange yet
//    for(i = 0; i < n; i++) {     // loop over rows to get scaling information
//       big = ZERO;
//       for(j = 0; j < n; j++) {
//          if((temp = fabs(a[i][j])) > big) big = temp;
//       }
//       if(big == ZERO) {
//          printf("\n\nSingular matrix in routine ludcmp()\n");
//          exit(1);
//       }               
//       vv[i] = 1.0/big;                 // save scaling */
//    } // end i-loop */

//    for(j = 0; j < n; j++) {     // loop over columns of Crout's method
//       for(i = 0; i< j; i++) {   // not i = j
//          sum = a[i][j];    
// 	 for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
// 	 a[i][j] = sum;
//       }
//       big = ZERO;   // initialization for search for largest pivot element
//       for(i = j; i< n; i++) {
//          sum = a[i][j];
// 	 for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
// 	 a[i][j] = sum;
// 	 if((dum = vv[i]*fabs(sum)) >= big) {
//   	    big = dum;
// 	    imax = i;
// 	 }
//       } // end i-loop
//       if(j != imax) {    // do we need to interchange rows ?
//          for(k = 0;k< n; k++) {       // yes
// 	    dum        = a[imax][k];
// 	    a[imax][k] = a[j][k];
// 	    a[j][k]    = dum;
// 	 }
// 	 (*d)    *= -1;            // and change the parit of d
// 	 vv[imax] = vv[j];         // also interchange scaling factor 
//       }
//       indx[j] = imax;
//       if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;

//         /*
//         ** if the pivot element is zero the matrix is singular
//         ** (at least to the precision of the algorithm). For 
//         ** some application of singular matrices, it is desirable
//         ** to substitute ZERO for zero,
//         */

//       if(j < (n - 1)) {                   // divide by pivot element 
//          dum = 1.0/a[j][j];
// 	 for(i=j+1;i < n; i++) a[i][j] *= dum;
//       }
//    } // end j-loop over columns
  
//    delete [] vv;   // release local memory

// }  // End: function ludcmp()




// void inverse(std::vector<std::vector<double>> a, int n){        
//   int          i,j, *indx;
//   double       d, *col, **y;

//   // allocate space in memory
//   indx = new int[n];
//   col  = new double[n];
//   std::vector<std::vector<double>> y(n, std::vector<double>(n, (double) 0)); 
   
//   ludcmp(a, n, indx, &d);   // LU decompose  a[][] 

//   // find inverse of a[][] by columns 

//   for(j = 0; j < n; j++) {

//     // initialize right-side of linear equations 

//     for(i = 0; i < n; i++) col[i] = 0.0;
//     col[j] = 1.0;

//     lubksb(a, n, indx, col);

//     // save result in y[][] 

//     for(i = 0; i < n; i++) y[i][j] = col[i];

//   }   //j-loop over columns 
   
//   // return the inverse matrix in a[][] 

//   for(i = 0; i < n; i++) {
//     for(j = 0; j < n; j++) a[i][j] = y[i][j];
//   } 
//   free_matrix((void **) y);     // release local memory 
//   delete [] col;
//   delete []indx;

// }  // End: function inverse()