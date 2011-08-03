/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BSEL WILL EVALUATE THE JDU, JDV-TH DERIVATIVE OF SURFACE B-REP. 
// ***INPUT* 
//    KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,NCD) 
//    ,ISR1,ISR2,NCD 
//      ....PROVIDE SURFACE B-REP, OF ORDER KU ALONG U AND KV ALONG V, 
//            U B-REP DIMENSION LUD, V B-REP DIMENSION LVD, 
//            AND SPACE DIMENSION NCD. 
//    U,V..........PARAMETER U AND V OF THE SURFACE B-REP. 
//    JDU,JDV......ARE ORDER OF DERIVATIVE ALONG U AND V DIRECTION. 
//              JDU,JDV MUST BE NON-NEGATIVE, MAY BE ZERO. 
// ***OUTPUT* 
//    P(NCD)    THE JDU,JDV-TH DERIVATIVE(S) EVALUATED . 
// *** NOTE * 
//    SURF(I,J,L),1<=I<=LUD 1<=J<=LVD, CONSTRUCT ONE B-COEF'S OF 
//    ONE COORDINATES. 1<=L 
void bse_(
	int ku, int lud, const double *ukt, int kv, int lvd, const double *vkt,
	const double *surf, int isr1, int isr2, int ncd,
	double u, double v, int jdu, int jdv, double *p);
