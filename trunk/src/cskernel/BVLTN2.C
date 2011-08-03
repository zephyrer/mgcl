/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Ble.h"
#include "cskernel/blgint.h"
#include "cskernel/bkdtkt.h"

// BVLATAN GETS APPROXIMATE VECTOR OF JDERIV-TH DERIVATIVE, 
// GIVEN POINT SEQUENCE P(.,.) WITH DATA POINTS TAU(.). 
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
//   JDERIV..ORDER OF DERIVATIVE. 
// ***OUTUT* 
//   DX(NCD)....VECTOR OF JDERIV-TH OBTAINED. 
//   IFLAG....=1: SUCCESSFUL RETURN 
//            =2: SOME ERROR DETECTED AND DX(.) WAS NOT CALCULATED 
// ***WORK* 
//   WORK(70),WKT(11),WRCOEF(8,NCD).... OF EACH LENGTH 
void bvltn2_(int np,const double *p,const double *tau, 
	int ip, int ncd, int ise, int jderiv,
	double *work, double *wkt, double *wrcoef, double *dx, int *iflag)
{
    // Local variables 
    int j;
    double s;
    int n2, istrt, kp;
    int kp2;

    // Parameter adjustments 
    --tau;
    --dx;
    p -= ip + 1;
    --work;

    // Function Body 
    for(j=1; j<=ncd; ++j)
		dx[j] = 0.f;
    *iflag = 2;
    if(np<=1)
		return;

    kp = np;
    if(kp>4)
		kp=4;
    kp2 = kp+2;
    if(kp2>np)
		kp2=np;
    if(ise==1){
		istrt = 1;
		s = tau[1];
    }else{
		istrt=np-kp2+1;
		s=tau[np];
    }
	// Here KP2 is the number of points to use of the point seq. P. 
	// Compute 1st derivative without end  condition. 
    n2 = kp2;
    bkdtkt_(&tau[istrt], n2, kp, wkt);
    *iflag=blgint_(&tau[istrt], &p[istrt + ip], wkt, kp, n2, ncd, ip,8, &work[kp2 + 1], wrcoef);
    if(*iflag != 1)
		return;
    ble_(kp, n2, wkt, wrcoef, 8, ncd, s, jderiv, &dx[1]);
}
