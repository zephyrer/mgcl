/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkktdp.h"
#include "cskernel/b1a.h"
#include "cskernel/blgint.h"
#include "cskernel/bk2fli.h"
#include "cskernel/bkmlt.h"
#include "cskernel/Ble.h"

//  REAL FUNCTION TO GET INTGERAL OF F(T)=X*DY-Y*DX. THIS INTEGRAL IS 
// FOR AREA COMPUTATION, I.E. AREA BOUNDED BY B-REP. 
// *** INPUT * 
// K,N,T(N+K),RCOEF(IRC,2).....B-REP OF ORDER K AND 2 SPACE DIMENSION. 
// T1,T2.....PARAMETER RANGE OF THE INTEGRATION, FROM T1 TO T2. 
// *** WORK * 
// TW(MM+2*K-2),WK1(MM),WK2(MM,4*K-4)....WORK AREA OF EACH LENGTH. HERE 
//         MM=N IF K=2, OTHERWISE MM=(N-K)*K+4*K-4. 
// *** OUTPUT * 
// BLA2D....THE INTEGRAL OF THE FUNCTION (X*DY-Y*DX), HERE X=RCOEF(.,1), 
//         Y=RCOEF(.,2), DX=DX/DT, DY=DY/DT. 
double bla2d_(int k, int n,const double *t,const double *rcoef, 
	int irc, double t1, double t2, double *tw, 
	double *wk1, double *wk2)
{
    // System generated locals 
    int rcoef_offset, isp1;
    double ret_val;

    // Local variables 
    int mltf, knew, nnew;
    double f[2];
    int i, l, iflag;
    int l1, l2;
    double df[2];
    int ie, is;
    int km1, km2;
    int npk, nep1, iem1;

// ****************** START OF BLA2D ***********************************
    // Parameter adjustments 
    --t;
    rcoef_offset = irc + 1;
    rcoef -= rcoef_offset;
    --tw;
    --wk1;
    --wk2;

    // Function Body 
    ret_val = 0.f;

    km2 = k - 2;
    km1 = k - 1;
    npk = n + k;
    l1=bk2fli_(npk, &t[1], t1);
    l2=bk2fli_(npk, &t[1], t2);
    if(l1 > l2){
		l = l1; l1 = l2;
		l2 = l;
    }
    ie = l1 - km2; if(ie<2) ie=2;
    nep1 = l2+1; if(nep1 > n+1) nep1= n+1;
    knew = (k << 1) - 2;

do{
	if(knew==2){
		is=ie; ie=nep1;
	}else{
		is = ie + km2;
	    isp1 = is + 1;
		bkmlt_(nep1, &t[1], isp1, km1, &ie, &mltf);
	    if(ie == -1) ie=nep1;
	}
// ===== 1. OBTAIN KNOT OF THE FUNCTION (X*DY/DT-Y*DX/DT).===== 
// ....   (1) FIND STARTING AND ENDING KNOT LOCATION. 
	while(t[is] == t[is+1]) is++;
	while(t[ie] == t[ie-1]) --ie;
    if(ie <= is) break;

    nnew = 0;
	iem1=ie-1;
// ..... (2) SET STARTING knew-KNOTS 
    for(l = 1; l <= knew; ++l) tw[++nnew] = t[is];
// ..... (3) SET INTERNAL KNOTS 
	i = is;
	if(knew==2){
		while(i<iem1) tw[++nnew] = t[++i];
	}else{
		while(i<iem1){
			do tw[++nnew] = t[++i]; while(t[i+1]==t[i]);
		    for(l = 1; l <= km1; ++l) tw[++nnew] = t[i];
		}
	}
// ..... (4) SET ENDING knew-KNOTS. 
    for(l=1; l<=knew; ++l) tw[++nnew] = t[ie];

    nnew -= knew;
//    KNEW=NEW ORDER, AND NNEW=NEW B-REP DIMENSION. 
// ===== 2. COMPUTE DATA POINTS AND THE ASSOCIATED DATA OF F(T). ===== 
    bkktdp_(nnew, knew, &tw[1], &wk1[1]);
    for(i = 1; i <= nnew; ++i){
		ble_(k,n,&t[1],&rcoef[rcoef_offset],irc,2,wk1[i],0,f);
		ble_(k,n,&t[1],&rcoef[rcoef_offset],irc,2,wk1[i],1,df);
		wk2[i] = f[0] * df[1] - f[1] * df[0];
    }// L400: 

//	for(i=1; i<=nnew; i++) printf("%d=%f,",i,wk1[i]);
//	for(i=1; i<=nnew+knew; i++) printf("%d=%f,",i,tw[i]);
// ===== 3. B-REP OF (X*DY-Y*DX). 
    iflag=blgint_(&wk1[1], &wk2[1], &tw[1], knew, nnew, 1, 1, 1, &wk2[nnew+1], &wk2[1]);
    if(iflag != 1) break;
// ===== 4. INTEGRAL OF (X*DY-Y*DX) ===== 
    ret_val += b1a_(knew, nnew, &tw[1], &wk2[1], t2) 
				- b1a_(knew, nnew, &tw[1], &wk2[1], t1);
}while(ie<nep1);

    return ret_val;
}
