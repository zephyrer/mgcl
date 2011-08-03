/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
void bvi2pl_(const double *g1,const double *g2, double *sl, int *iflag);
