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
double bvs3pd_(const double *v1,const double *v2,const double *v3);
