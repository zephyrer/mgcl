/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// INTERNAL SUBROUTINE OF BLGCS, 
// BLGCS7 REMOVES DATA POINTS THAT ARE TOO CLOSE TO PREVIOUS DATA POINT 
// *** INPUT * 
//     TAU(N)...DATA POINT OF LENGTH N 
//     VAL(IV,NCD)..DATA POINT ORDINATE OF LENGTH N, 
//     KWORK(N).....INDICATES (TAU(I),VAL(I,.)) MAY BE REMOVED OR NOT: 
//          KWORK(I)=1 NOT ALLOWED TO REMOVE,  =0 ALLOWED. 
//     N........NUM OF DATA POINTS. 
//     NCD......SPACE DIMENSION OF THE DATA VAL 
//     IV.......ROW DIMENSION OF VAL 
// *** OUTPUT * 
//     TAU(N)....DATA POINTS UPDATED 
//     VAL(IV,NCD)...DATA POINT ORDINATES UPDATED 
//     N.....NEW NUM OF DATA POINTS 
void blgcs7_(double *tau, double *val, int *kwork, int *n, int ncd, int iv);
