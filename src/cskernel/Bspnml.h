/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BSPNML GETS NORMAL VECTOR OF SURFACE B-REP. 
// *** INPUT  * 
//     KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,3) 
//            ...SURFACE B-REP 
//     U,V.....PARAMETER VALUES OF THE ABOVE SURFACE B-REP. 
// *** OUTPUT * 
//     VNML(3) : NORMAL VECTOR 
void bspnml_(int ku, int lud,const double *ukt, 
	int kv, int lvd,const double *vkt,const double *surf, int isr1, int isr2,
	double u, double v, double *vnml);
