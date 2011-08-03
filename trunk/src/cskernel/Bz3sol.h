/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BZ3SOL solves linear equations of 3 variables, i.e. 
// solves the system  A*X=B, where A and B are known and X is unknown. 
// *** INPUT *** 
//   A(3,3)........Left hand side 3 by 3 matrix. 
//   B(3)..........Right hand side 3 by 1 matrix. 
// *** OUTPUT *** 
//   X(3)........solution. 
void bz3sol_(const double *a,const double *b, double *x,int *iflag);