/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVVRT computes the vector Y, such that three vector E, Y, and E*X 
// construct right hand orthogonal system. This means the vector Y is 
// perpendicular to the plane which includes both vectors E and X. 
// *** INPUT *** 
//   E(3),X(3)..... are two input vectors. 
// *** OUTPUT *** 
//   Y(3).......... is output vector. 
// *** NOTE *** 
// Two vectors E and X must not be linear, this causes divide abend. 
void bvvrt_(const double *e,const double *x, double *y);