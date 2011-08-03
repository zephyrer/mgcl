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
void bvvpd_(const double *x,const double *y, double *z);