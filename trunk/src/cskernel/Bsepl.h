/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
	double *rcoef);
