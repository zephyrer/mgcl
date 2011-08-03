/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLGCS2 is a dedicated subroutine of BLGCS, generates data point 
// sequence TAU(.) from data (VAL(j,NCD),j=1...N). VAL(j,.) may 
// include circle inf that is declared by KVAL(j). 
// ******Input****** 
//   NCD,N,KVAL(N),VAL(IV,NCD).....Input data of Space Dimension NCD, and 
//         of length N. KVAL(j) is a knuckle inf of VAL(j,.). 
// ******Output***** 
//   TAU(N).........Data point abssisa obtained of length N. 
void blgcs2_(int ncd, int n, const int *kval, 
	const double *val, int iv, double *tau
);
