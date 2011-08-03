/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLG4S1 IS A DEDICATED SUBROUTINE OF BLG4SQ. 
// BLG4S1 PRODUCES SHOENBERG'S VARIATION DIMINISHING KNOTS OF ORDER 4. 
// DATA POINTS TAU MAY BE MULTIPLE, I.E. 
// WHEN     TAU(I-1) < TAU(I)=---=TAU(I+M-1) < TAU(I+M)  , 
//          T(I+2+J) = TAU(I+J) FOR J=0,1,--,M-1 ( M<=3) . 
// *** INPUT * 
//     k..........order of the b-spline.
//     TAU(N).....DATA POINT SEQUENCE. 
//     N..........NUMBER OF DATA POINTS. 
// *** OUTPUT* 
//     T(N+4).....KNOT VECTOR OF SHOENBERG'S VARIATION DIMINISHING. 
//     IFLAG......INDICATES SUCCESS (=1) , OR FAILURE (=2). 
//          FAILURE IS CAUSED BY ILLEGAL NUMBER OF D.P. MULTIPLICITY; 
//           1) MORE THAN 3 MULTIPLICITY DETECTED. 
//           2) NO MUTIPLICITY ON STARTING (ENDING) AND NEXT (PREVIOUS) 
//              DATA POINT IS MULTIPLE. 
void blg4s1_(int k, const double *tau, int n, double *t, int *iflag);
