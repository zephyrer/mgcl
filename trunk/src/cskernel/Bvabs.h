/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVABS computes length of a vector VECTOR(.,NCD) 
// INPUT *** 
//            NCD ..... is space dimension of the VECTOR. 
//            IVEC..... is the row dimension of VECTOR, i.e. 
//                     VECTOR is declared as VECTOR(IVEC,NCD) 
//            VECTOR... the vector to evaluate. 
// OUTPUT *** 
//            BVABS....vector length of VECTOR. 
double bvabs_(int ncd, int ivec,const double *vector);
