/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blg4s1.h"
#include "cskernel/blg4sp2.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//  mgblgsq IS EASY-TO-USE VERSION OF blg4sp2_, I.E. mgblgsq GENERATES 
//  KNOT VECTOR USING blg4s1_, THEN B-COEF'S USING blg4sp2_. 
// *** INPUT * 
//         k......order of the b-spline.
//         IBCBEG, IBCEND......BOUNDARY COND OF BEGINNING AND ENDING 
//                POINT, EACH. 
//                =1 1ST DERIV PROVIDED,    =2 2ND DERIV 
//                =3 NO BOUNDARY COND. 
//                =4 BOTH 1ST AND 2ND DDERIVATIVES PROVIDED. 
//         TAU(N) : DATA-POINTS 
//         VAL(IV,NCD) : KOSIN-OFFSET-DATA OF M-ORDER 
//         IV     : COLUMN-LENGTH OF VAL 
//         N      : DATA-NO. OF VAL      ( IV >= N ) 
//         NCD    : ORDER OF COORDINATES 
//         IRC    : COLUMN-LENGTH OF RCOEF 
// *** OUTPUT * 
//         T(N+k) : KNOT VECTOR 
//         RCOEF(IRC,NCD) : B-SPLINE 
//         IFLAG  : =1 NORMAL END 
//                  =2 ABNORMAL 
// *** WORK * 
//         WORK(N,k*2+1) : WORK AREA FOR SUBROUTINE blg4sp2_ 

// ****    CREATE KNOT-VECTOR T(N+4) FROM TAU(N) 
int mgblgsq(int k, int ibcbeg, int ibcend, double *tau,
	double *val, int iv, int n, int ncd, int irc,
	double *work, double *t, double *rcoef, int* iflag)
{
	int i;
    // Parameter adjustments 
    work -= n+1; val -= iv+1; rcoef -= irc+1;

    // Function Body 
    blg4s1_(k,tau, n, t, iflag);
    if (*iflag == 1) {
	// ****    CREATE B-SPLINE RCOEF(N,NCD) FROM TAU(N),T(N+4),VAL(IV,NCD)  
		*iflag = 2;
		for (i = 1; i <= ncd; ++i) {
		    blg4sp2_(k,iflag, ibcbeg, ibcend, tau, &val[i*iv+1], 
			    iv, n, 1, t, 1, &work[n+1], &work[(n<<1) + 1], 
				&work[n*3+1], &rcoef[i*irc+1]);
			if (*iflag == 2) return 0;
		}
    }
    return 0;
}
