/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/bvvpd.h"
#include "cskernel/bvabs.h"

// FUNCTION TO OBTAIN DISTANCE BETWEEN a LINE(G) AND a POINT(P). 
// INPUT *** 
//    NCD......SPACE DIMENSION OF G AND P, MUST BE 2 OR 3. 
//    G(NCD,2) : PARAMETER OF S.L. AS BELOW 
//            G(.,1):A POINT ON THE S.L., G(.,2):DIRECTIONAL COSINE 
//    P(NCD) : COORDINATE OF A POINT 
// OUTPUT *** BVDPSL : DISTANCE 
double bvdpsl_(int ncd,const double *g,const double *p){
    double ret_val;
    double a[3], c[3], e[3];
    int j;
    double r;

    // Parameter adjustments 
    --p;
    g -= ncd+1;;

    // Function Body 
    a[2] = 0.f;
    e[2] = 0.f;
    for(j=1; j<=ncd; ++j){
		a[j-1] = p[j]-g[j+ncd];
		e[j-1] = g[j+(ncd<<1)];
    }
    bvvpd_(a,e,c);
    bvvpd_(e,c,c);
    r = bvabs_(3,1,c);
    if(r<=bzmzro_())
		ret_val = 0.f;
    else
		ret_val = (c[0]*a[0]+c[1]*a[1]+c[2]*a[2])/r;
    return ret_val;
}
