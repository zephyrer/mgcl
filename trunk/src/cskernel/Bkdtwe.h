/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//     BKDTWE WILL GENERATE DATA POINT SEQUENCE IN TAU(I), 1 <=I<=N, 
//     S.T.          TAU(1) = 0. 
//                   TAU(I) = DIST(P(I,.)-P(I-1,.)) + TAU(I-1) 
//                             FOR 2<=I<=N . 
//     IF IBCBEG,IBCEND <> 3 , MULTIPLE D.P. ADDED AT THE ASSOCIATED 
//     D.P. 
// *** INPUT * 
//     IBCBEG,IBCEND INDICATES BOUNDAY POINT INF ,I.E. 
//              =3:ALL DATA OF VAL ARE POSITIONAL DATA 
//             <>3:VAL(2,.) OR VAL(N-1,.) ARE DERIVATIVE INF, AND 
//                 MULTIPLE D.P. ADDED AT THE POINT 
//     VAL(IV,NCD)   POINT-SEQUENCE IN A NCD-DIMENSION SPACE. 
//     N             NUM OF POINTS, P(I,.) 1<=I<=N  ARE EFFECTIVE. 
//     NCD           SPACE DIMENSION OF POINTS P( , ). 
//                   ( COLUMN DIMENSION OF THE VARIABLE P. ) 
//     IV            ROW DIMENSION OF THE VARIABLE VAL 
// *** OUTPUT * 
//     TAU(N)        DATA POINT SEQUENCE. 
void bkdtwe_(int ibcbeg, int ibcend,const double *val, int n, int ncd, int iv, double *tau);
