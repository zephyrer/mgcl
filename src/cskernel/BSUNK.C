/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Blunk.h"

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
){
    int i, j;
    int jj, kk;

    // Parameter adjustments 
    surf1 -= is11*(is12+1)+1;
    work2 -= lvd2+1;
    surf2 -= is21*(is22+1)+1;

	// --- U DIR. AND BOTH DIR. 
    if(iknot==1 || iknot==3){
		for(kk=1; kk<=3; ++kk){
		for(j=1; j<=lvd1; ++j)
			blunk_(ku,lud1,ukt1,&surf1[(j+kk*is12)*is11+1],
				lvd2,1,lud2,ukt2,lvd2,work1,&surf2[(j+kk*is22)*is21+1]);
		}
		// --- V DIR. 
		if(iknot==3){
			for(kk=1; kk<=3; ++kk){
			for(i=1; i<=lud2; ++i){
		    for(j=1; j<=lvd1; ++j)
				work2[j+lvd2]=surf2[i+(j+kk*is22)*is21];
				blunk_(kv,lvd1,vkt1,&work2[lvd2+1],
					lvd2,1,lvd2,vkt2,lvd2,work1,&work2[(lvd2<<1)+1]);
		    for(jj=1; jj<=lvd2;++jj)
				surf2[i+(jj+kk*is22)*is21]=work2[jj+(lvd2<<1)];
			}
			}
		}

    }else if(iknot==2){
		for(kk=1; kk<=3; ++kk){
		for(i=1; i<=lud1; ++i){
		for(j=1; j<=lvd1; ++j)
		    work2[j+lvd2]=surf1[i+(j+kk*is12)*is11];
		blunk_(kv,lvd1,vkt1,&work2[lvd2+1],
			lvd1,1,lvd2,vkt2,lvd2,work1,&work2[(lvd2<<1)+1]);
		for(jj=1; jj<=lvd2; ++jj)
		    surf2[i+(jj+kk*is22)*is21]=work2[jj+(lvd2<<1)];
		}
		}
    }
}
