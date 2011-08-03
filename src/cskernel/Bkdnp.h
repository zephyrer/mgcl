/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BKDNP deletes too near data point(s). 
//
// ***** Input ****** 
//  N,TAU(N), VAL(IV,NCD) 
//        Data points of length N, 
//        (TAU(i),VAL(i,j)) is the data point sequence, must be 
//        non-decreasing. TAU(i) is 
//        the abssisa, and VAL(i,j) is the ordinates, 1 <= i <= N. 
//        NCD is the space dimension, and 1 <= j <= NCD. 
//        IV is the row dimension of the variable VAL. 
// IMLT   Indicates whether multiple data points be deleted or not. 
//        Multple data points mean the points, 
//           TAU(i)=TAU(i+1). 
//        IMLT =1  : delete multiple points. 
//             <>1 : not delete multiple data points. 
// RATIO  Ratio to regard as a too near point. ratio>1.
//        Let r = (TAU(i)-TAU(i-1))/(TAU(i-1)-TAU(i-2)), then 
//        if r > RATIO, TAU(i-1) is removed. 
//        if r < 1/RATIO, TAU(i-1) is removed.
// Let dtold=tau[i-1]-tau[i-2], and dtnew=tau[i]-tau[i-1], then
// if dnew>ratio*dtold, tau[i-2] will be removed.
// if dtold>dtnew*ratio, tau[i] will be removed.
//
// ***** Output ***** 
// N,TAU(N), VAL(IV,NCD)   will be updated. 
//
// ***** NOTE ***** 
// Start and End Points are never removed, instead the next point is 
// removed. 
void bkdnp_(int *n, double *tau, double *val, int iv, int ncd, int imlt, double ratio);
