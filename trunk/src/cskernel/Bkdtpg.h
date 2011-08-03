/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//     bkdtpg_ WILL GENERATE DATA POINT SEQUENCE IN tau[I], 0 <=I<=n-1, 
//     S.T.          tau[0] = 0. 
//                   tau[I] = DIST(p(I,.)-p(I-1,.)) + tau[I-1] 
//                             FOR 1<=I<=n-1 . 
// *** INPUT * 
//     p(ip,ncd)     POINT-SEQUENCE IN A ncd-DIMENSION SPACE. 
//     n             NUM OF POINTS, p(I,.) 1<=I<=n  ARE EFFECTIVE. 
//     ncd           SPACE DIMENSION OF POINTS p( , ). 
//                   ( COLUMN DIMENSION OF THE VARIABLE p. ) 
//     ip            ROW DIMENSION OF THE VARIABLE p . 
// *** OUTPUT * 
//     tau(n)        DATA POINT SEQUENCE. 
void bkdtpg_(const double *p, int n, int ncd, int ip, double *tau);
