/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/bkdtpg.h"
#include "cskernel/Bvunit.h"
#include "cskernel/Bvltan.h"
#include "cskernel/bkdnp.h"

// BVUTAN GETS APPROXIMATE UNIT TANGENT VECTOR, GIVEN POINT SEQUENCE. 
// THE VECTOR MAY BE STARTING AT P(1,.), OR ENDING AT P(NP,.) 
// ***INPUT* 
//   ISE....INDICATES WHETHER STARTING OR ENDING VECTOR 
//           =1: STARTING AT P(1,.)    =2: ENDING AT P(NP,.) 
//   NP....GIVES NUM OF INPUT POINTS 
//   P(IP,NCD)....PROVIDES POINT SEQUENCE OF NCD SPACE DIMENSION 
//   IP....IS ROW DIMENSION OF THE VARIABLE P 
//   NCD...IS SPACE DIMENSION OF TNPUT POINTS 
// ***OUTUT* 
//   RTAN(NCD).....THE UNIT TANGENT VECTOR 
//   IFLAG....=1: SUCCESSFUL RETURN 
//            =2: SOME ERROR DETECTED AND RTAN(.) WAS NOT CALCULATED 
// ***WORK* 
//   WORK(81+NCD*8) 
void bvutan_(int ise, int np,const double *p, 
	int ip, int ncd, double *work, double *rtan, 
	int *iflag)
{
    int i, j, n;
    int istrt;
    double tau[6], pnt[18];	//was [6][3]

    // Parameter adjustments 
    p -= ip+1;

    // Function Body 
    n=np;
    if(n>6)
		n = 6;
   
    istrt=0;
    if(ise==2)
		istrt = np - n;
    for(i=1; i<=n; ++i){
		for(j=1; j<=ncd; ++j)
			pnt[i+j*6-7] = p[istrt+i+j*ip];
    }
    bkdtpg_(pnt, n, ncd, 6, tau);
    bkdnp_(&n, tau, pnt, 6, ncd, 1, bkmax_());

	// GET TANGENT VECTOR 
    bvltan_(n,pnt,tau,6,ncd,ise,work,work+70,work+81,rtan,iflag);
    bvunit_(rtan, ncd, rtan);
}
