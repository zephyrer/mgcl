/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//         BLIP2 = MAX(ABS(DATA(MU+I)-F)) FOR I=0,1,---,NDATA-1 
// BLIP2 is for BLIPP3. 
// *** INPUT * 
//     RDATA(MU+NDATA-1)  : DATA SEQUENCE 
//     MU                : THE SUBSCRIPT OF DATA WHICH INDICATES WHERE 
//                         BLIP2 SHOULD START. 
//     NDATA             : NUMBER OF DATA TO HANDLE 
// *** OUTPUT * 
//     BLIP2 = MAX(ABS(RDATA(MU+I)-F)) FOR I=0,1,---,NDATA-1 
double blip2_(const double *rdata, int mu, int ndata, double f);
