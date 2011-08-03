/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// bkcrng_ WILL CONVERT THE OLD KNOT SEQUENCE t[] INTO tnew[] SO THAT 
// tnew[k-1]=tstrt, tnew[n]=tend, AND MULTIPLICITY OF tnew[I] AND tnew[n+I] 
// HOLD. (0<=I<=k-1). 
//   WHEN k=0, t IS ASSUMED AS DATA POINT SEQUENCE. 
// *** INPUT * 
//     n,k,t[n+k]....INPUT KNOT VECTOR OF A B-REP, 
//                   n:B-REP DIMENSION, k:ORDER, t[.]:KNOT VECTOR 
//                   k MAY BE 0, AND INDICATES t[.] ARE DATA POINTS. 
//     tstrt,tend.......NEW PARAMETER RANGE OF NEW KNOT VECTOR,I.E. 
//                   VALUES SUCH THAT tnew[k-1]=tstrt AND tnew[n]=tend 
// *** OUTPUT * 
//     tnew[n+k].....NEW KNOT VECTOR OBTAINED 
//                   tnew MAY BE THE SAME AREA AS t. 
void bkcrng_(int k, int n, const double *t, double tstrt, double tend, double *tnew);
