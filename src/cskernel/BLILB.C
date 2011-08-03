/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blipp.h"

int blilb_(int *k, int *n, double *t, double 
	*rcoef, int *irc, double *g, double *error, double *
	ts, int *mx, double *work1, double *rcoefy, int *nx, 
	double *x, int *iend)
{
    /* System generated locals */
    int rcoef_dim1, rcoef_offset, i__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    double d__, c__, f;
    int i__;
    double s;

/* BLILB WILL GET PARAMETER VALUES X(I),I=1..NX ,OF INTERSECTION POINTS */
/* OF A STRAIGHT LINE G(.) AND B-REP.  X(.)>TS AND NX<=MX. */
/* *** INPUT  * */
/*     K,N,T(N+K),RCOEF(IRC,2),IRC.....B-REP OF 2 SPACE DIMENSION. */
/*     G(2,2).......VECTOR TO EXPRESS THE STRAIGHT LINE */
/*              G(.,1) :  A POINT ON THE LINE */
/*              G(.,2) :  THE DIRECTIONAL VECTOR */
/*                  THE LINE IS EXPRESSED USING PARAMETER T AS, */
/*                  X=G(1,1)+G(1,2)*T, Y=G(2,1)+G(2,2)*T */
/*     ERROR        ALLOWABLE ERROR VALUE TO GET INTERSECTION POINTS */
/*     TS.......START PARAM VALUE OF THE INTERSECTION COMPUTATION. */
/*     MX.......INDICATES LENGTH OF THE VARIABLE X(.), I.E. X(MX). */
/* *** OUTPUT * */
/*     NX           NUMBER OF THE SOLUTIONS. IF NX=O, NO PARAM. FOUND. */
/*     X(NX)        PARAMETER VALUES OBTAINED BY BLILB */
/*     IEND.........INDICATES WHETHER COMPUTATION REACHED TO THE END. */
/*                  =0: NOT TO THE END,   =1:TO THE END. */
/* *** WORK   * */
/*     WORK1(4*K*K+3K),RCOEFY(N) WORK ARRAY OF EACH LENGTH */

/* ****** START OF BLILB ****** */
/*       X'=X*COS(THETA)+Y*SIN(THETA) */
/*       Y'=-X*SIN(THETA)+Y*COS(THETA) */

    /* Parameter adjustments */
    --rcoefy;
    --t;
    rcoef_dim1 = *irc;
    rcoef_offset = rcoef_dim1 + 1;
    rcoef -= rcoef_offset;
    g -= 3;
    --x;
    --work1;

    /* Function Body */
    d__ = sqrt(g[5] * g[5] + g[6] * g[6]);
    c__ = g[5] / d__;
    s = g[6] / d__;
/*     CHANGE OF RCOEF TO RCOEFY */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	rcoefy[i__] = -rcoef[i__ + rcoef_dim1] * s + rcoef[i__ + (rcoef_dim1 
		<< 1)] * c__;
    }
/*     CALCULATION OF F (FUNCTION VALUE TO GET THE INTERSECTION POINTS) */
/*     CHANGE OF THE STRAIGHT LINE BY ROTATION */
    f = -g[3] * s + g[4] * c__;
/*     GET PARAMETER VALUE OF B-REP. USING BLIPP. */
    blipp_(*k,*n,&t[1],&rcoefy[1],f,*error,*ts,*mx, &work1[1], nx, &x[1], 
	    iend);
    return 0;
} /* blilb_ */

