/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkktdp.h"
#include "cskernel/bspnml.h"
#include "cskernel/blgint.h"

// BSTPER COMPUTES TANGENT PLANE OF ONE OF GIVEN SURFACE B-REP PERIMETER. 
// *** INPUT *** 
// KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,3)...SURFACE B-REP .
//  KR......PERIMETER NUM:    1     2     3     4 
//                          V=MIN  U=MAX V=MAX U=MIN 
// *** OUTPUT *** 
//  K,N,T(N+K),RCOEF(IRC,3).....T.P. B-REP OBTAINED (ORDER K=MAX(KU,KV)) 
// *** WORK ARRAY *** 
//  WORK(N,K*2): WORK AREA 
void bstper_(int ku, int lud,const double *ukt, 
	int kv, int lvd,const double *vkt,const double *surf, int isr1, int isr2, int ncd,
	int kr, int irc, double *work, int *k, int *n, double *t, double *rcoef)
{
    int ircp1;
    int i, j, iflag, npk;
    double vec[3], tsd;

    // Parameter adjustments 
    ircp1 = irc + 1;
    rcoef -= ircp1;

    // Function Body 
    if(kr==1 || kr==3){
		*n = lud;
		*k = ku;
		bkktdp_(*n, ku, ukt, work);// *** DATA POINT GEN. *** 
		if(kr==1)
			tsd=vkt[kv-1];
		else
		    tsd=vkt[lvd];
		// *** NORMAL VECTOR GEN. *** 
		for(i=1; i<=*n; ++i){
		    bspnml_(ku,lud,ukt,kv,lvd,vkt,surf,isr1,isr2,work[i-1],tsd,vec);
			for(j=1; j<=3; ++j)
				rcoef[i+j*irc] = vec[j-1];
		}
		//      *** KNOT GEN. *** 
		npk = *n + *k;
		for(i=0; i<npk; ++i)
			t[i] = ukt[i];
    }else{
		*n = lvd;
		*k = kv;
		bkktdp_(*n,kv,vkt,work);// *** DATA POINT GEN. *** 
		if(kr==4)
			tsd = ukt[ku-1];
		else
		    tsd = ukt[lud];
		//      *** NORMAL VECTOR GEN. *** 
		for(i=1; i<=*n; ++i){
		    bspnml_(ku,lud,ukt,kv,lvd,vkt,surf,isr1,isr2,tsd,work[i-1],vec);
			for(j=1; j<=3; ++j)
				rcoef[i+j*irc] = vec[j-1];
		}
		//      *** KNOT GEN. *** 
		npk = *n + *k;
		for(i=0; i<npk; ++i)
			t[i] = vkt[i];
    }
    iflag=blgint_(work,&rcoef[ircp1],t,*k,*n,3,irc,irc,work+*n,&rcoef[ircp1]);
}
