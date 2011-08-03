/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVCSL computes the center of a circle that osculates to two input 
// straight lines. 
// *** INPUT *** 
// R..... Radius of the circle. 
// NCD... is the space dimension, 2 or 3. 
// P(NCD),G1(NCD),G2(NCD).......are two stright lines that pass 
// point P(.) and their directional vectors are G1(.) and G2(.). 
// *** OUTPUT *** 
// CENTER(NCD).....center of the circle 
// P1(NCD),P2(NCD),T1,T2..... Are the osculating points at the straight 
// lines, P1 is the point coordinate of 1st straight line (P,G1) and T1 
// is the parameter value of the straight line expression. (P2,T2) are 
// the same as to 2nd straight line. 
void bvcsl_(double r, int ncd,const double *p, 
	const double *g1,const double *g2,double *center, double *p1, 
	double *p2, double *t1, double *t2);
