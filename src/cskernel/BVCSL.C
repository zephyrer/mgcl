/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bvvpd.h"
#include "cskernel/Bvunit.h"
#include "cskernel/bvinpd.h"
#include "cskernel/bvi2sl.h"

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
	double *p2, double *t1, double *t2)
{
    int j;
    double a1[3], a2[3], q1[3], q2[3];
    double c2d[2];
    double rnorml[3], g12d[2], g22d[2], p12d[2], p22d[2], pns1[3], pns2[3];

// <<<<< 1. Get straight lines tha are moved normal to the originals.>>>>> 
// (P,G1) AND (P,G2) re the original line and (PNS1,G1), (PNS2,G2) 
// are the new lines obtained. 

    // Parameter adjustments 
    --p2;
    --p1;
    --center;
    --g2;
    --g1;
    --p;

    // Function Body 
    q1[2] = 0.f;
    q2[2] = 0.f;
    for(j=1; j<=ncd; ++j){
		q1[j-1] = g1[j];
		q2[j-1] = g2[j];
    }
    bvvpd_(q1,q2,rnorml);
    bvvpd_(rnorml,q1,a1);
    bvvpd_(q2, rnorml,a2);
    bvunit_(a1,ncd,a1);
    bvunit_(a2,ncd,a2);
    for(j=1; j<=ncd; ++j){
		pns1[j-1] = p[j]+r*a1[j-1];
		pns2[j-1] = p[j]+r*a2[j-1];
    }

//<<<<< 2. Get the center of the osculating circle and contact points.>>>>>
    if(ncd==2){
		bvi2sl_(pns1,&g1[1],pns2,&g2[1],&center[1],t1,t2);
    }else{
		bvunit_(q1, 3, q1);
		p12d[0] = bvinpd_(3,pns1,q1);
		p12d[1] = bvinpd_(3,pns1,a1);
		g12d[0] = bvinpd_(3,&g1[1],q1);
		g12d[1] = 0.f;
		p22d[0] = bvinpd_(3,pns2,q1);
		p22d[1] = bvinpd_(3,pns2,a1);
		g22d[0] = bvinpd_(3,&g2[1],q1);
		g22d[1] = bvinpd_(3,&g2[1],a1);
		bvi2sl_(p12d,g12d,p22d,g22d,c2d,t1,t2);
		for(j=1; j<=3; ++j)
		    center[j] = c2d[0]*q1[j-1]+c2d[1]*a1[j-1];
    }
    for(j=1; j<=ncd; ++j){
		p1[j] = p[j]+g1[j]*(*t1);
		p2[j] = p[j]+g2[j]*(*t2);
    }
}
