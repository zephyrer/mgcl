/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include <malloc.h>
#include "cskernel/tolerance.h"
#include "cskernel/bpval2.h"
#include "cskernel/blcbp.h"

#define MAX_ORDER 20
#define MAX_LOOP 50

// REAL FUNCTION TO GET INTERSECTION POINT OF 1-DIMENSIONAL B-REP. 
// B-REP. IS FIRST CONVERTED INTO PP-REP, THEN NEWTON-RAPHSON METHOD IS
//  USED FOR THE SOLUTION. EVALUATION OF B-REP IS DONE USING PP-REP. 
// *** INPUT *
//   K,N,T(N+K),RCOEF(N).......B-REP FOR INTERSECTION COMPUTATION. 
//   KI........KNOT INDEX OF T S.T. 
//                        BLER(,T(KI),) <= F <= BLER(,T(KI+1),)   , 
//                    OR  BLER(,T(KI),) >= F >= BLER(,T(KI+1),) 
//   F.......THE B-VALUE TO FIND THE ASSOCIATED PARAMETER VALUE BLI1SP, 
//                                  F = BLER(,BLI1SP,). 
//         ERROR               ERROR ESTIMATE OF F, I.E. THE VALUE S.T. 
//  ERROR......ERROR ALLOWED TO COMPUTE INTERSECTION. 
// *** OUTPUT * 
//    BLI1SP........THE PARAMETER VALUE S.T.   F=BLER(,BLI1SP,) 
// *** WORK * 
//   WORK(K*K)         WORK AREA OF K*K 
double bli1sp_(int k, int n, const double *t, const double *rcoef, 
	int ki, double f, double error, double *work){
    int i,is,nloop;
    double d1,tdif,rbrk[2],fl, fm, fr;
    double tl, tm, tr, d2l, d2r, zero, tzero;
    double pcoef_buf[MAX_ORDER];

// **************************** START OF BLI1SP ************************
	double* pcoef=pcoef_buf;
	if(k>MAX_ORDER)
		pcoef=(double*)malloc(sizeof(double)*k);

	if(ki>n)
		ki=n;
    while(t[ki-1] >= t[ki])
		++ki;
    
// CONVERT B-REP INTO PP-REP. 
    is = ki-k;
    blcbp_(k,k,&t[is],&rcoef[is],work,rbrk,pcoef,&i);
//  NOW PP-REP IS ACHIEVED. 
    tl = rbrk[0];
    tr = rbrk[1];
	tzero=bzrzro_()*(tr-tl);

	fl = pcoef[0]-f;
    fr = bpval2_(rbrk,pcoef,k,tr,1,0)-f;

	if(fl*fr>=0.){
		if(fl>=fr){
			if(fl<=0.){tm=tl; goto L8888;}
			if(fr>=0.){tm=tr; goto L8888;}
		}else{
			if(fl>=0.){tm=tl; goto L8888;}
			if(fr<=0.){tm=tr; goto L8888;}
		}
	}

	//  NOW  FL < F < FR(or FL > F > FR), GET THE SOLUTION BY NEWTON-RAPHSON METHOD. 
	zero = bzmzro_();
	while(1){
		if(tr-tl<=tzero)	goto L8888;
		if(fabs(fl)<=error){tm=tl; goto L8888;}
		if(fabs(fr)<=error){tm=tr; goto L8888;}
		d2l = bpval2_(rbrk,pcoef,k,tl,1,2);
		d2r = bpval2_(rbrk,pcoef,k,tr,1,2);
		if(d2l*d2r>=0.){
			if(d2r*fr>0.){
				tm = tr;
				fm = fr;
			}else{
			    tm = tl;
				fm = fl;
			}
			// GET 1ST DERIV AT tm FOR ITERATION. 
		    d1 = bpval2_(rbrk,pcoef,k,tm,1,1);
			if(fabs(d1)>zero){
				tdif = fm/d1;
				tm -= tdif;
				if(tl<=tm && tm<=tr){

			// ITERATE UNTIL ABS(F-FM) < ERROR . 
				for(nloop=0; nloop<MAX_LOOP; nloop++){
					if (fabs(tdif) <= fabs(error/d1)) goto L8888;
					fm = bpval2_(rbrk,pcoef,k,tm,1,0)-f;
					if(fabs(fm)<=error) goto L8888;
					d1 = bpval2_(rbrk,pcoef,k,tm,1,1);
					tdif = fm/d1;
					tm -= tdif;
					if(tm<tl || tr<tm) break;
					if(fabs(tdif)<=tzero) goto L8888;
				}
			
				}
			}
		}

		// NEWTON-RAPHSON METHOD WILL NOT CONVERGE, NARROW THE INTERVAL BY 
		// BI-SECTION METHOD. 
		tm = (tl+tr)*0.5;
		fm = bpval2_(rbrk,pcoef,k,tm,1,0)-f;
		if(fl*fm<0.){
			fr=fm; tr=tm;
		}else{
			tl=tm; fl=fm;
		}

	}

L8888:
	if(k>MAX_ORDER) free(pcoef);
    return tm;
}
