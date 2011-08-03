/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bvvpd.h"
#include "cskernel/bvinpd.h"
#include "cskernel/bvcang.h"

// BVIVEC interpoates two vectors A and B, given interpolation 
// ratio R. 
// *** INPUT *** 
//   NCD..... is the space dimension of the two vector A and B. 
//   A(NCD), B(NCD)......two vectors to interpoate. 
//   R........is the ratio between 0 and 1, indicating how much A 
//           be rotated. When R=0, X=A, and when R=1, X=B. 
// *** OUTPUT *** 
//   X(NCD)...The unit vector interpolating cector A and B. 
//           Ratio of THETA1 and THETA2 is R, where THETA1 and THETA2 
//           are the angles between A nad X, and B and X, each. 
//   T1,T2...Parameter value to express X, as 
//           X=T1*A+T2*B. 
// *** NOTE* 
//   Vector A and B mut be unit vectors. 
//   NCD must be 2 or3. 
// Table of constant values 
void bvivec_(int ncd,const double *a,const double *b, double r,
			double *x, double *t1, double *t2){
    static double thtmin = .5;

    int i,j;
    double theta, theta1, theta2;
    double xx, aa, bb;
    double e1[3], e2[3], e3[3], vzero;

    // Parameter adjustments 
    --x;
    --b;
    --a;

    vzero = bzmzro_();

	// Function Body 
    aa = 0.f;
    bb = 0.f;
    for (j = 1; j <= ncd; ++j) {
		aa += a[j]*a[j];
		bb += b[j]*b[j];
    }
    aa = sqrt(aa);
    bb = sqrt(bb);

    if (aa <= vzero || bb <= vzero) {
		*t2 = r;
		*t1 = 1.f - r;
    }else{
		e1[2]=e2[2]=0.f;//Initialize for the case of ncd==2.
		for (j = 1; j <= ncd; ++j) {
		    e1[j-1] = a[j]/aa;
			e2[j-1] = b[j]/bb;
		}
		bvvpd_(e1, e2, e3);
		xx = sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
 		if (xx < thtmin) {
		    theta = asin(xx);
			if (bvinpd_(ncd,&a[1],&b[1]) < 0.f)
				theta = 3.1415926535897932384626433833f - theta;
		} else {
			xx = bvcang_(ncd,&a[1],&b[1]);
		    theta = acos(xx);
		    xx = sin(theta);
		}
		theta1 = r*theta;
		theta2 = theta-theta1;
		*t1 = sin(theta2)/xx;
		*t2 = sin(theta1)/xx;
	}

    for (i = 1; i <= ncd; ++i) {
		x[i] = a[i] * *t1 + b[i] * *t2;
    }
}
