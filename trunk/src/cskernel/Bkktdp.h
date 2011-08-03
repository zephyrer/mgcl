/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BKKTDP WILL GENERATE DATA POINT SEQUENCE TAU(I), 0<=I<=N-1, FROM 
// KNOT VECTOR OF A B-REP. 
// *** INPUT * 
//     N,K,T(N+K)....INPUT KNOT VECTOR OF A B-REP, 
//                   N:B-REP DIMENSION, K:ORDER, T(.):KNOT VECTOR 
// *** OUTPUT * 
//     TAU(N).....DATA POINT SEQ. GENERATED. 
void bkktdp_(int n, int k, const double *t, double *tau);
