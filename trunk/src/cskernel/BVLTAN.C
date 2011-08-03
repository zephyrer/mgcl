/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Bvltn2.h"

// BVLATAN GETS APPROXIMATE TANGENT VECTOR, GIVEN POINT SEQUENCE P(.,.) 
// WITH DATA POINTS TAU(.). 
// THE VECTOR MAY BE STARTING AT P(1,.), OR ENDING AT P(NP,.) 
// ACCORDING TO THE FLAG, ISE. 
// ***INPUT* 
//   NP....GIVES NUM OF INPUT POINTS 
//   P(IP,NCD)....PROVIDES POINT SEQUENCE OF NCD SPACE DIMENSION 
//   TAU(NP)....DATA POINT ABSISSA OF P(.) 
//   IP....IS ROW DIMENSION OF THE VARIABLE P 
//   NCD...IS SPACE DIMENSION OF TNPUT POINTS 
//   ISE....INDICATES WHETHER STARTING OR ENDING VECTOR 
//           =1: STARTING AT P(1,.)    <>1: ENDING AT P(NP,.) 
// ***OUTUT* 
//   DX(NCD)....TANGENT VECTOR OBTAINED. 
//   IFLAG....=1: SUCCESSFUL RETURN 
//            =2: SOME ERROR DETECTED AND DX(.) WAS NOT CALCULATED 
// ***WORK* 
//   WORK(70),WKT(11),WRCOEF(8,NCD).... OF EACH LENGTH 
void bvltan_(int np,const double *p,const double *tau, 
	int ip, int ncd, int ise, double *work, double *wkt, double *wrcoef, double *dx, int *iflag
){
	bvltn2_(np,p,tau,ip,ncd,ise,1,work,wkt,wrcoef,dx,iflag);
}
