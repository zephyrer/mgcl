/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blipp.h"

void blibp_(int k, int n, double *t, double	*rcoef, int irc,
			double *g, double error, double ts, int mx, double *work, double *rcoefx,
			int *nx, double *x, int *iend)
{
    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    double d;
    int i;
    double gl[4];

/*     BLIBP WILL GET ALL PARAMETER VALUES X(I),I=1,NX, OF INTERSECTION */
/*     POINTS OF A PLANE WITH B-REP. */
/* *** INPUT  * */
/*   K,N,T(N+K),RCOEF(IRC,3),IRC....PROVIDE 3-D SPACE DIMENSION B-REP. */
/*      G(4)...........EXPRESS THE PLANE AS BELOW: */
/*                     G(1)*X+G(2)*Y+G(3)*Z=G(4) */
/*      ERROR          ALLOWABLE ERROR VALUE TO GET INTERSECTION POINTS */
/*      TS.............START PARAM VALUE OF INTERSECTION COMPUTATION, */
/*                     I.E. X(I) > TS FOR 1<=I<=NX<MX. */
/*      MX.............SIZE OF VARIABLE X(.). */
/* *** OUTPUT * */
/*      NX             NUMBER OF THE SOLUTIONS */
/*      X(NX)          PARAMETER VALUES WHICH ARE OBTAINED BY BLIBP */
/*      IEND...........INDICATES WETHER BLIBP COMPUTED TO THE END OF */
/*                     B-REP(IEND=1), OR NOT TO THE END(IEND=0) BECAUSE */
/*                     MX IS TOO SMALL. */
/* *** WORK   * */
/*      WORK(4*K*K+3*K),RCOEFX(N)......WORK ARRAY OF EACH LENGTH. */
/*  ***START OF BLIBP *** */
    /* Parameter adjustments */
    rcoef -= irc + 1;
    --rcoefx;

    /* Function Body */
    d = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
    for(i=0; i<4; ++i)
		gl[i] = g[i]/d;

/*     CHANGE OF CO-ORDINATE SYSTEM BY ROTATION */
/*      X'=X*G(1)+Y*G(2)+Z*G(3) */
/*     THEN,PLANE WILL BECOME PARARELL TO X'-AXIS. */
    for(i=1; i<=n; ++i)
		rcoefx[i] = rcoef[i+irc]*gl[0]+rcoef[i+(irc<<1)]*gl[1]+rcoef[i+irc*3]*gl[2];
/*     FUNCTION VALUE TO GET THE INTERSECTION POINTS */
/*     GET PARAMETER VALUES OF B-REP. USING BLIPP */
    blipp_(k,n,t,&rcoefx[1],gl[3],error,ts,mx,work,nx,x,iend);
}
