/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bvvpd.h"
#include "cskernel/Bvunit.h"
#include "cskernel/bse.h"

// BSPNML GETS NORMAL VECTOR OF SURFACE B-REP. 
// *** INPUT  * 
//     KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,3) 
//            ...SURFACE B-REP 
//     U,V.....PARAMETER VALUES OF THE ABOVE SURFACE B-REP. 
// *** OUTPUT * 
//     VNML(3) : NORMAL VECTOR 
void bspnml_(int ku, int lud,const double *ukt, 
	int kv, int lvd,const double *vkt,const double *surf, int isr1, int isr2,
	double u, double v, double *vnml)
{
    double p01[3], p10[3];

    // Function Body 
	bse_(ku,lud,ukt,kv,lvd,vkt,surf,isr1,isr2,3,u,v,1,0,p10);
    bse_(ku,lud,ukt,kv,lvd,vkt,surf,isr1,isr2,3,u,v,0,1,p01);
    bvvpd_(p10,p01,vnml);
    bvunit_(vnml,3,vnml);
}
