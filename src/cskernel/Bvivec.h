/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
			double *x, double *t1, double *t2);
