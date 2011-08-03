/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// SUBROUTIN TO OBTAIN SUBDIVIDED B-COEFFICIENTS WITH NEW KNOT 
// INPUT*** IKNOT  : KIND OF SUBDIVIDED KNOT 
//                   =1 SUBDIVIDE ONLY UKT1( U-DIR. ) 
//                   =2 SUBDIVIDE ONLY VKT1( V-DIR. ) 
//                   =3 SUBDIVIDE BOTH UKT1 AND VKT1 ( BOTH U,V-DIR. ) 
//          KU     : ORDER ALONG U-DIRECTION 
//          LUD1   : OLD DIMENSION OF B-REP. ( U DIR. ) 
//          UKT1(LUD1+KU) : OLD KNOT VECTOR OF SURF. B-REP. IN U DIR. 
//          KV     : ORDER ALONG V-DIRECTION 
//          LVD1   : OLD DIMENSION OF B-REP. ( V DIR. ) 
//          VKT1(LVD1+KV) : OLD KNOT VECTOR OF SURF. B-REP. IN V DIR. 
//          SURF1(IS11,IS12,3) : OLD B-COEFFICIENT OF SURF B-REP. 
//          LUD2   : NEW DIMENSION OF B-REP. ( U DIR. ) 
//          UKT2(LUD2+KU) : NEW KNOT VECTOR OF SURF. B-REP. IN U DIR. 
//          LVD2    : NEW DIMENSION OF B-REP. ( V DIR. ) 
//          VKT2(LVD2+KV) : NEW KNOT VECTOR OF SURF. B-REP. IN V DIR. 
// OUTPUT *** SURF2(IS21,IS22,3) : NEW B-COEFFICIENT OF SURF B-REP. 
// WORK ***   WORK1(K,K),WORK2(LVD2,2) WHERE K=MAX(KU,KV) 
void bsunk_(int iknot, int ku, int lud1,const double *ukt1,
	int kv, int lvd1,const double *vkt1, 
	const double *surf1, int is11, int is12,
	int lud2,const double *ukt2, int lvd2,const double *vkt2,
	int is21, int is22, double *work1, double *work2, double *surf2
);
