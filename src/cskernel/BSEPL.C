/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Bleval.h"

// BSEPL COMPUTES PARAMETER LINE, GIVEN U OR V PARAMETER VALUE. 
// *** INPUT * 
//     KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,NCD) 
//      ,ISR1,ISR2,NCD 
//         ....SURFACE B-REP OF ORDER KU,KV AND NCD SPACE DIMENSION. 
//     KX,X....INDICATE WHICH PARAMETER LINE BE OBTAINED AS: 
//              KX=1  :X=U-VALUE,   V-PARAMETER LINE 
//              KX<>1 :X=V-VALUE,   U-PARAMETER LINE 
//     JDERIV..SPECIFIES ORDER OF DERIVATIVE OF THE PARAMETER LINE. 
//     IRC.....ROW DIMENSION OF RCOEF 
// *** OUTPUT * 
//     K,N,T(N+K),RCOEF(IRC,NCD)..LINE B-REP OBTAINED.  N=LUD OR LVD 
//             ACCORDING TO KX. THE ORDER K IS KU OR KV. 
void bsepl_(int ku, int lud, const double *ukt, 
	int kv, int lvd, const double *vkt, const double *surf, int isr1, int isr2, int ncd,
	int kx, double x, int jderiv, int irc, int *k, int *n, double *t, 
	double *rcoef)
{
    int i, j, npk;

    surf -= isr1*(isr2+1)+1;;
    rcoef -= irc+1;
    if(kx==1){
	// ***** V - PARAMETER LINE 
		*k = kv;
		npk = lvd+kv;
		for(i=0; i<npk; ++i)
			t[i] = vkt[i];
		*n = lvd;
		for(j=1; j<=ncd; ++j){
			bleval_(ku,lud,ukt,&surf[(j*isr2+1)*isr1+1],isr1,
				lvd,kx,x,jderiv,1,j,&rcoef[j*irc+1]);
		}
    }else{
	// ***** U - PARAMETER LINE 
		*k = ku;
		npk = lud + ku;
		for(i=0; i<npk; ++i)
		    t[i] = ukt[i];
		*n = lud;
		for(i=1; i<=lud; ++i){
	    for(j=1; j<=ncd; ++j){
			bleval_(kv,lvd,vkt, &surf[i+(j*isr2+1)*isr1],isr1,
				1,kx,x,jderiv,1,j,&rcoef[i+j*irc]);
		}
		}
    }
}
