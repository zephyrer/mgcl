/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"

// BVI2SL computes the intersection point of two 2D straight lines, 
// (P1(2),DC1(2)) and (P2(2),DC2(2)). 
// **** INPUT *** 
// P1(2),DC1(2)......1st straight line, P1(2) is a point and DC1(2) is 
//          the directional cosine. 
//          The straight line can be expressed as follows using 
//          parameter t; X(t)=P1(1)+DC1(1)*t 
//                       Y(t)=P1(2)+DC1(2)*t 
// P2(2),DC2(2)......is 2nd straight line. 
// **** OUTPUT *** 
//    P(2),T1,T2......are intersection point coordinate data P(.) 
//          and the parameter value of 1st(T1) and 2nd line(T2). 
void bvi2sl_(const double *p1,const double *dc1,const double *p2,const double *dc2,
			double *p, double *t1, double *t2){
    double r, x1, x2;
    double zerror;

    // Parameter adjustments 
    --dc2;
    --p2;
    --dc1;
    --p1;

    // Function Body 
    zerror = bzamin_();
    r = dc1[2]*dc2[1]-dc1[1]*dc2[2];
    if(fabs(r)<=zerror){
		*t1 = 0.f;
		*t2 = 0.f;
    }else{
		x1 = p2[1]-p1[1];
		x2 = p2[2]-p1[2];
		*t1 = (-dc2[2]*x1+dc2[1]*x2)/r;
		*t2 = (-dc1[2]*x1+dc1[1]*x2)/r;
    }
    p[0] = p1[1]+dc1[1]*(*t1);
    p[1] = p1[2]+dc1[2]*(*t1);
}
