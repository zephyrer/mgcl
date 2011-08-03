/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/b2dv.h"
#include "cskernel/bk2fli.h"

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
void bse_(int ku, int lud, const double *ukt, int kv, int lvd, const double *vkt,
	const double *surf, int isr1, int isr2, int ncd, double u, double v, int jdu, int jdv, double *p)
{
    // Local variables 
    int lefu, lefv, lumk, lvmk, i__, j, m;
    double rcoef;
    int ludpk, lvdpk;
    int ii, jj;

// **************** START OF BSEL ***********************************
    double rbatu_buf[20], rbatv_buf[20];
	double* rbatu=rbatu_buf; double* rbatv=rbatv_buf;
	if(ku>20)
		rbatu=(double*)malloc(sizeof(double)*ku);
	if(kv>20)
		rbatv=(double*)malloc(sizeof(double)*kv);

    // Parameter adjustments 
    --p;
    surf -= isr1*(isr2+1)+1;

    ludpk = lud+ku;
    lvdpk = lvd+kv;
//      GET COEFFICIENTS OF B-COEF'S. 
//       U-DIRECTION 
    lefu=bk2fli_(ludpk, ukt, u);
    b2dv_(ku, ludpk, ukt, u, lefu, jdu, rbatu);
//       V-DIRECTION 
    lefv=bk2fli_(lvdpk, vkt, v);
    b2dv_(kv, lvdpk, vkt, v, lefv, jdv, rbatv);

    lumk = lefu - ku;
    lvmk = lefv - kv;
//   NOW COEFFICIENTS OF B-COEF'S ARE OBTAINED IN RBATU ,RBATV. 
//   EVALUATION 
    for(m=1; m<=ncd; ++m){
		p[m] = 0.f;
		for(j=1; j<=kv; ++j){
			rcoef = 0.f;
		    jj = lvmk + j;
		    if (jj < 1)
				jj = 1;
		    if (jj > lvd)
				jj = lvd;
			for(i__ = 1; i__ <= ku; ++i__){
				ii = lumk + i__;
				if(ii<1)
				    ii = 1;
				if (ii > lud)
				    ii = lud;
				rcoef += rbatu[i__-1]*surf[ii+(jj+m*isr2)*isr1];
		    }
			p[m] += rbatv[j-1]*rcoef;
		}
    }
	if(ku>20) free(rbatu);
	if(kv>20) free(rbatv);
}
