/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/bvabs.h"

// BVCANG EVALUATES COSINE ANGLE OF TWO VECTORS. LET THETA IS 
// THE ANGLE OF TWO VECTORS, THEN BVCANG=COS(THETA). 
// INPUT *** 
//     NCD...... is the space dimension of the two vector V1 
//               and V2. 
//     IV1,IV2....ROW DIMENSION OF VECTOR V1 AND V2. 
//     V1(IV1,NCD),V2(IV2,NCD)........ are two vectors, input as 
//               VECTOR V1(1,.) and VECTOR V2(1,.). 
// OUTPUT *** BVCANG : COSINE ANGLE OF THE TWO VECTOR. 
double bvcang_(int ncd,const double *v1,const double *v2){
    double ret_val;
    int i;
    double r;
    double r1, r2, vzero;
    double v1v2;

    r1 = bvabs_(ncd,1,v1);
    r2 = bvabs_(ncd,1,v2);
    r = r1*r2;
    vzero = bzmzro_();
    vzero *= vzero;
    if(r<=vzero)
		ret_val = 0.f;
    else{
		v1v2 = 0.f;
		for(i=0; i<ncd; ++i)
			v1v2 += v1[i]*v2[i];
		ret_val =v1v2/r;
    }
    return ret_val;
}
