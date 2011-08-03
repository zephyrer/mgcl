/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bk1fli.h"

//bpval_ CALCULATES VALUE AT X OF JDERIV-TH DERIVATIVE OF
//PP Function FROM PP-REPR.
//
// ******  I N P U T  ****** 
//  tau, coef, L, K.....FORMS THE PP-REPRESENTATION OF THE FUNCTION  F 
//        TO BE EVALUATED. SPECIFICALLY, THE J-TH DERIVATIVE OF  F  IS 
//        GIVEN BY 
//     (D**J)F(X) =coef(J+1,I) + H*(coef(J+2,I) + H*( ... (coef(K-1,I) + 
//                             + H*coef(K,I)/(K-J-1))/(K-J-2) ... )/2)/1 
//        WITH  H = X - tau(I),  AND 
//       I  =  MAX( 1 , MAX( J ,  tau(J) .LE. X , 1 .LE. J .LE. L ) ). 
//  X.....THE POINT AT WHICH TO EVALUATE. 
//  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUAT- 
//        ED.  A S S U M E D  TO BE ZERO OR POSITIVE. 
// ******  O U T P U T  ****** 
//  BPVAL.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF  F  AT  X. 
// ******  M E T H O D  ****** 
//     THE BKINTRAL INDEX  I , APPROPRIATE FOR  X , IS FOUND THROUGH A 
//  CALL TO  BK1FLI . THE FORMULA ABOVE FOR THE  JDERIV-TH DERIVATIVE 
//  OF  F  IS THEN EVALUATED (BY NESTED MULTIPLICATION). 
double bpval_(const double *tau,const double *coef, int l, int k, double x, int jderiv){
    int i, m;
    double fmmjdr, ret_val, h;

	if(k<=jderiv)
		return 0.;
		//DERIVATIVES OF ORDER  K  OR HIGHER ARE IDENTICALLY ZERO. 

    // Parameter adjustments 
    --tau;
    coef -= k+1;

//FIND INDEX  I  OF LARGEST BREAKPOINT TO THE LEFT OF  X . 
    i=bk1fli_(l, &tau[1], x);
    if(i<1)
		i = 1;

    ret_val = 0.f;
    fmmjdr = (double) (k-jderiv);
//EVALUATE  JDERIV-TH DERIVATIVE OF  I-TH POLYNOMIAL PIECE AT  X . 
    h = x-tau[i];
	for(m=k; m>jderiv; m--){
	    ret_val = ret_val/fmmjdr*h + coef[m+i*k];
	    fmmjdr += -1.f;
	}
    return ret_val;
}
