/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/Bvunit.h"
#include "cskernel/blipp.h"
#include "cskernel/bk1fli.h"

//                                                     Y. MIZUNO 
// BLUMI2 IS AN INTERNAL SUBROUTINE OF BLUMIX, 
// BLUMI2 WILL GET PARAM VALUE OF B-REP (N,T,RCOEF) WHOSE FUNCTION VALUE 
// IS F(.). 
// *** INPUT  * 
//     KCOD1,K,N,T(.),RCOEF(IRC,2),IRC....DESCRIBE THE B-REP OF ORDER K. 
//     KCOD2,TAU,F(2)......ARE DATA POINT AND ASSOCIATED FUNCTION VALUE 
//            , KCOD2 SPECIFIES COORDINATE KIND. 
//     TPREV.....IS THE PARAM VALUE USED WHEN MORE THAN TWO INTERSECTION 
//            POINTS, NEAREST AND GREATER THAN TPREV PARAM VALUE IS 
//            EMPLOYED. 
// *** OUTPUT *** 
//     TINT......PARAM VALUE OF THE B-REP OBTAINED. 
// ***WORK * WORK(4*K*K+3K)....WORK OF LENGTH 4*K*K+3K. 
void blumi2_(int kcod1, int k, int n, 
	const double *t, const double *rcoef, int irc, int kcod2, 
	double tau, double *f, double tprev, double error, 
	double *work, double *tint)
{
    int rcoef_offset;

    // Local variables 
    int iend, left, kcom1, kcom2;
    double x[10];
    int id, nx;
    double ts;

    // Parameter adjustments 
    rcoef_offset = irc + 1;
    rcoef -= rcoef_offset;
    --t;
    --f;

    // Function Body 
    ts = t[k]-1.f;
    if(kcod1==4)
	//   CASE OF GIRTH B-REP 
		blipp_(k,n,&t[1],&rcoef[rcoef_offset],tau,error,ts,10,work,&nx,x,&iend);
    else{
		//   CASE OF NOT GIRTH B-REP 
		if (kcod2 == 4) {
		//      CASE OF F(1) IS GIRTH 
			*tint = f[1];
		    return;
		}else{
		    kcom1 = 1;
		    if(kcod1%3+1==kcod2)
				kcom1=2;
	  
			kcom2=kcom1%2+1;
		    blipp_(k,n,&t[1],&rcoef[kcom1*irc+1],f[kcom2],error,ts,10,work,&nx,x,&iend);
		    if (nx == 0) {
				blipp_(k,n,&t[1],&rcoef[kcom1*irc+1],f[kcom2]-error,error,ts,10,work,&nx,x,&iend);
				if(nx==0)
				    blipp_(k,n,&t[1],&rcoef[kcom1*irc+1],f[kcom2]+error,error,ts,10,work,&nx,x,&iend);
			}
		}
    }

    if(nx==0){
		*tint=t[n+1];
    }else if(nx<=1){
		*tint=x[0];
    }else{
		left=bk1fli_(nx,x,tprev);
		id = left+1;
		if(id>nx)
			id = nx;
		*tint = x[id-1];
    }
}

// BLUMI4         BLUMI4 FOR BLUMIX 
// BLUMI4 GETS UNIT TANGENT VECTOR FROM TWO TWO-DIMENSIONAL TANGENT 
// VECTORS F1(2) AND F2(2) AND STORE THEM IN RCOEF(1,J) FOR 1<=J<=3. 
//   *** KCOD1 MUST NOT BE 4 *** 
void blumi4_(int kcod1, double *f1, int kcod2, 
	double *f2, double *rcoef, int irc)
{
    double vect[3];
    int kcom1, kcom2, j, kcod11, kcod12;

    // Parameter adjustments 
    --f1;
    --f2;
    rcoef -= irc + 1;

    // Function Body 
    kcod11 = kcod1%3+1;
    kcod12 = kcod11%3+1;
    vect[kcod11-1] = f1[1];
    vect[kcod12-1] = f1[2];

    if(kcod2!=4){
		kcom1 = 1;
		if(kcod11==kcod2)
		    kcom1 = 2;
		kcom2 = kcom1 % 2 + 1;
		if(fabs(f2[kcom2]) > bzmzro_()) {
			vect[kcod1 - 1] = f2[kcom1] * (f1[kcom1] / f2[kcom2]);
		} else {
		    vect[kcod1 - 1] = f2[kcom1];
		}
	// CASE OF GIRTH B-REP 
    }else{
		vect[kcod1-1] = f2[2];
		vect[kcod11-1] = f1[1] * f2[1];
		vect[kcod12-1] = f1[2] * f2[1];
    }
	//   CONVERT INTO NIT VECTOR 
    bvunit_(vect,3,vect);
    for(j=1; j<=3; ++j)
		rcoef[j*irc+1] = vect[j-1]*1.008f;
}
