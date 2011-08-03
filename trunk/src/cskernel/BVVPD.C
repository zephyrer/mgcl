/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVVPD computes vector product Z of two vector X and Y, 
// i.e. Z=X*Y. 
// *** INPUT *** 
//   X(3), Y(3)..... are two input 3D vectors. 
// *** OUTPUT *** 
//   Z(3).......... is output vector. 
void bvvpd_(const double *x,const double *y, double *z){
    z[0] = x[1]*y[2] - x[2]*y[1];
    z[1] = x[2]*y[0] - x[0]*y[2];
    z[2] = x[0]*y[1] - x[1]*y[0];
}
