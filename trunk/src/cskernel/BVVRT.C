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
void bvvrt_(const double *e,const double *x, double *y){
    // Builtin functions 
    double sqrt(double);
    int i;
    double ee, xe, yy;

    ee = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    xe = x[0]*e[0] + x[1]*e[1] + x[2]*e[2];
    for(i=0; i<3; ++i)
		y[i] = ee*x[i] - xe*e[i];

    yy = y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
    yy = sqrt(yy);
    for(i=0; i<3; ++i)
		y[i] /= yy;
}
