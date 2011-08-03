/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// FUNCTION TO EVALUATE SCALOR TRIPLE PRODUCT OF VECTOR V1,V2,V3. 
//bvs3pd_=determinant(V1,V2,V3).
// ****INPUT**** 
//     V1(3),V2(3),V3(3).....THREE VECTORS OF SPACE DIMENSION 3. 
// ****OUTPUT*** 
//     BVS3PD.................TRIPLE PRODUCT(SCOLAR VALUE). 
double bvs3pd_(const double *v1,const double *v2,const double *v3){
    // Local variables 
    double w1, w2, w3;
    w1 = v1[0]*v2[1]*v3[2] - v1[2]*v2[1]*v3[0];
    w2 = v1[2]*v2[0]*v3[1] - v1[0]*v2[2]*v3[1];
    w3 = v1[1]*v2[2]*v3[0] - v1[1]*v2[0]*v3[2];
    return  w1+w2+w3;
}
