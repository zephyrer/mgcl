/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVINPD computes inner product of two vectors V1 and V2. 
// *** INPUT *** 
//   NCD,V1(NCD),V2(NCD)..... are two input vectors. 
//         NCD IS THE SPACE DIENSION. 
// *** OUTPUT *** 
//   BVINPD........is scalar data of the inner product. 
double bvinpd_(int ncd,const double *v1,const double *v2);
