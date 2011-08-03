/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// bkdtkt_ WILL GENERATE KNOT VECTORt from the data points tau.
// t[i], 0<=i<=n+k-1 are generated as:
//    t[0]=t[1]=--=t[k-1]=tau[0] 
//    t[i+k]=(tau[i+1]+---+tau[i+k-1])/(k-1) for 0<=i<=n-k-1
//    t[n]=t[n+1]=---=t[n+k-1] = tau[n-1] 
// *** INPUT * 
//     tau[n]    :  DATA POINT SEQUENCE 
//     n         :  NUM OF DATA POINTS. 
//     k         :  ORDER OF B-REP. 
// *** OUTPUT * 
//     t[n+k]    :  KNOT VECTOR 
void bkdtkt_(const double *tau, int n, int k, double *t);
