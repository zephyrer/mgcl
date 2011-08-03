/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// TO BE CALLED IN  B L G S M T 
// FROM A PRACTICAL GUIDE TO SPLINE BY C.DE-BOOR, AND 
// UPDATED BY Y. MIZUNO. 
//
// CONSTRUCTS THE UPPER THREE DIAGS. IN V(I,J), I=2,,NPOINT-1, J=1,3, OF 
//  THE MATRIX  6*(1-P)*Q-TRANSP.*(D**2)*Q + P*R, THEN COMPUTES ITS 
//  L*L-TRANSP. DECOMPOSITION AND STORES IT ALSO IN V, THEN APPLIES 
//  FORWARD AND BACK SUBSTITUTION TO THE RIGHT SIDE Q-TRANSP.*Y IN  QTY 
//  TO OBTAIN THE SOLUTION IN  U . 
// *** INPUT *** 
// P..................PARAMETER P. 
// V(NPOINT,7)........OUTPUT OF BLGSM1. 
// QTY(NPOINT,NCD)....Q-TRANSP. * Y (OUTPUT OF BLGSM1) 
// NPOINT.............NUM OF POINTS. 
// DY(I).........ERROR ESTIMATE, 1<=I<=NPOINT. 
// NCD................SPACE DIMENSION 
// *** OUTPUT *** 
// U(NPOINT,NCD) 
// QU(IQU,NCD)........Q*U 
// DQU................DY*QU. 
void blgsm2_(double p, double *v, const double *qty, 
	int npoint,const double *dy, int ncd, int iqu, 
	double *u, double *qu, double *dqu)
{
    // Local variables 
    double prev, dquw, dqux, twop;
    int i, j;
    double ratio, six1mp;
    int nm1, nm2, n2, n3, n4;

    // Parameter adjustments 
    --dy;
    v -= npoint+1;
    u -= npoint+1;
    qty -= npoint+1;
    qu -=  iqu + 1;;
    nm1 = npoint - 1;
    nm2 = npoint - 2;
	n2=npoint*2; n3=npoint*3; n4=npoint*4;
	six1mp = (1.f - p) * 6.f;
    twop = p * 2.f;

//     CONSTRUCT 6*(1-P)*Q-TRANSP.*(D**2)*Q  +  P*R 
    if(npoint<=2){
	    for(j=1; j<=ncd; ++j){
			u[j*npoint+1] = 0.f;
			u[j*npoint+2] = 0.f;
	    }
	}else{
	    for(i=2; i<=nm1; ++i){
			v[i+npoint] = six1mp*v[i+npoint*5]+twop*(v[i-1+n4]+v[i+n4]);
			v[i+n2] = six1mp*v[i+npoint*6]+p*v[i+n4];
			v[i+n3] = six1mp*v[i+npoint*7];
		}
	    if(nm2<2) {
		    for(j=1; j<=ncd; ++j){
				u[j*npoint+1] = 0.f;
				u[j*npoint+2] = qty[j*npoint+2]/v[npoint+2];
				u[j*npoint+3] = 0.f;
			}
		}else{
			//  FACTORIZATION 
		    for(i=2; i<=nm2; ++i){
				ratio = v[i+n2]/v[i+npoint];
				v[i+1+npoint] -= ratio*v[i+n2];
				v[i+1+n2] -= ratio*v[i+n3];
				v[i+n2] = ratio;
				ratio = v[i+n3]/v[i+npoint];
				v[i+2+npoint] -= ratio*v[i+n3];
				v[i+n3] = ratio;
			}
			v[n3+1] = 0.f;

		    for(j=1; j<=ncd; ++j){
				//  FORWARD SUBSTITUTION 
				u[j*npoint+1] = 0.f;
				u[j*npoint+2] = qty[j*npoint+2];
				for(i=2; i<=nm2; ++i)
				    u[i+1+j*npoint] = qty[i+1+j*npoint]-v[i+n2]*u[i+j*npoint]-v[i-1+n3]*u[i-1+j*npoint];

				//  BACK SUBSTITUTION 
				u[npoint+j*npoint] = 0.f;
				u[nm1+j*npoint] /= v[nm1+npoint];
				for(i=nm2; i>1; --i)
					u[i+j*npoint] = u[i+j*npoint]/v[i+npoint]-u[i+1+j*npoint]*v[i+n2]-u[i+2+j*npoint]*v[i+n3];
			}
		}
	}

//  CONSTRUCT Q*U 
    *dqu = 0.f;
    for(j=1; j<=ncd; ++j){
		prev = 0.f;
		for(i=2; i<=npoint; ++i){
			qu[i+j*iqu]=(u[i+j*npoint]-u[i-1+j*npoint])/v[i-1+n4];
			qu[i-1+j*iqu] = qu[i+j*iqu]-prev;
		    prev = qu[i+j*iqu];
		}
		qu[npoint+j*iqu] = -qu[npoint+j*iqu];

		dquw = 0.f;
		for(i=1; i<=npoint; ++i){
			dqux = dy[i]*qu[i+j*iqu];
		    dquw += dqux*dqux;
		}
		*dqu += dquw;
    }
}
