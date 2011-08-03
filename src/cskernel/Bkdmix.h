/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BKDMIX WILL MIX DATA SEQ. T1 AND T2, STORED IN T. 
// *** INPUT * 
//  N1          NUM OF DATA SEQ T1. 
//  T1(N1)      DATA SEQ 1. 
//  N2          NUM OF DATA SEQ T2. 
//  T2(N2)      DATA SEQ 2. 
// *** OUTPUT * 
//  N           NUM OF NEW DATA SEQ T. 
//  T(N)        NEW DATA SEQUENCE. 
// *** NOTE * 
//    DATA POINT MULTIPLICITY IS ALLOWED ONLY IN T1. 
void bkdmix_(int n1, const double *t1, int n2, const double *t2, unsigned *n, double *t);
