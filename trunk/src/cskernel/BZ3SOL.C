/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bvs3pd.h"

//BZ3SO1 COMPUTES I-J COFACTOR OF THE MATRIX A(3,3)
double bz3so1_(const double *a, int i, int j){
    int k1, k2, l1, l2;
    a -= 4;

    k1 = i%3+1;
    k2 = k1%3+1;
    l1 = j%3+1;
    l2 = l1%3+1;
    return a[k1+l1*3]*a[k2+l2*3]-a[k1+l2*3]*a[k2+l1*3];
}

// BZ3SOL solves linear equations of 3 variables, i.e. 
// solves the system  A*X=B, where A and B are known and X is unknown. 
// *** INPUT *** 
//   A(3,3)........Left hand side 3 by 3 matrix. 
//   B(3)..........Right hand side 3 by 1 matrix. 
// *** OUTPUT *** 
//   X(3)........solution. 
void bz3sol_(const double *a,const double *b, double *x,int *iflag){
    double c[9];	// was [3][3]
    int i, j;
    double detmin;
    double det;

    // Parameter adjustments 
    --x;
    --b;

    // Function Body 
    detmin = bzmzro_();
    det = bvs3pd_(a,a+3,a+6);
	if(fabs(det)>detmin){
		for(i=1; i<=3; ++i){
		    for(j=1; j<=3; ++j)
				c[i+j*3-4] = bz3so1_(a,j,i)/det;
		}
		for(i=1; i<=3; ++i)
			x[i]= c[i-1]*b[1]+c[i+2]*b[2]+c[i+5]*b[3];
		*iflag = 1;
    }else
		*iflag = 2;
}
