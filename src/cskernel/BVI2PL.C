/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bvvpd.h"

// Subroutine to compute the straight line that is the intersection 
// of two plane G1 and G2. 
// INPUT *** 
//    G1(4),G2(4) ....describes two planes as 
//                   Gi(1)*X+Gi(2)*Y+Gi(3)*Z=Gi(4)  for i=1,2 
// OUTPUT *** 
//    SL(3,2)..the straight line that is the intersection of G1 and G2. 
//         The straight line can be expressed as below using a parameter 
//            t. 
//             X=SL(1,1)+SL(1,2)*t 
//             Y=SL(2,1)+SL(2,2)*t 
//             Z=SL(3,1)+SL(3,2)*t. 
// Note1.  Three vectors G1(.), G2(.), and SL(.,1) constitutes right hand
//           orthogonal system. 
//      Note2.  The vector SL(.,1) is a unit vector. 
//    IFLAG .........  Return code, 
//            =1: intersection line obtained normaly, 
//            =2: two planes are parallel. 
void bvi2pl_(const double *g1,const double *g2, double *sl, int *iflag){
    int j;
    double r, slabs[3];
    int i1, i2, i3;

    // Parameter adjustments 
    sl -= 4;
    --g2;
    --g1;

    // Function Body 
    bvvpd_(&g1[1], &g2[1], &sl[7]);
    r = sqrt(sl[7]*sl[7]+sl[8]*sl[8]+sl[9]*sl[9]);
    if(r>bzmzro_()){
		for(j=1; j<=3; ++j){
		    sl[j+6]/=r;
			slabs[j-1]=fabs(sl[j+6]);
		}
		i1 = 3;
		if(slabs[0]>=slabs[1]){
		    if(slabs[0]>=slabs[2])
				i1 = 1;
		}else{
		    if(slabs[1]>=slabs[2])
				i1 = 2;
		}
		i2 = i1%3+1;
		i3 = i2%3+1;
		r = g1[i3]*g2[i2]-g2[i3]*g1[i2];
		sl[i1+3] = 0.f;
		sl[i2+3] = (g2[4]*g1[i3]-g1[4]*g2[i3])/r;
		sl[i3+3] = (g1[4]*g2[i2]-g2[4]*g1[i2])/r;
		*iflag = 1;
    }else
		*iflag = 2;
}
