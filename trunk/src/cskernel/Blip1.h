/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLIP1               BLIP1 ESPECIALLY FOR BLIP 
//      BLIP1 WILL GET THE INTERSECTION NUMBER BETWEEN THE POLYGON 
//      POLYGN(I) AND THE STRAIGHT LINE F. 
// *** INPUT * 
//         POLYGN(MU+NSPAN)  : SPECIFIES POLYGON. 
//         MU                : INDICATES AT WHICH POINT OF POLYGN BLIP1 
//                             SHOULD START THE COUNTING. 
//         NSPAN             : INDICATES AT WHICH POINT OF POLYGN BLIP1 
//                             SHOULD STOP THE COUNTING. 
//         F                 : SPECIFIES THE STRAIGHT LINE. 
// *** OUTPUT * 
//         BLIP1            : NUMBER OF INTERSECTION POINT. 
int blip1_(const double *polygn, int mu, int nspan, double f);
