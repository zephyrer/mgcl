/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//TO BE CALLED IN B L G S M T 
//FROM A PRACTICAL GUIDE TO SPLINE AND UPDATED BY Y.MIZUNO,9/7,'83 
//
//PUT DELX = X(.+1) - X(.)  INTO  V(.,4), 
//PUT THE THREE BANDS OF  Q-TRANSP*D  INTO  V(.,1-3), AND 
//PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q)  AT AND ABOVE THE DIAGONAL 
//INTO  V(.,5-7) . 
//HERE,  Q IS  THE TRIDIAGONAL MATRIX OF ORDER (NPOINT-2,NPOINT) 
//WITH GENERAL ROW  1/DELX(I) , -1/DELX(I) - 1/DELX(I+1) , 1/DELX(I+1) 
//AND   D  IS THE DIAGONAL MATRIX  WITH GENERAL ROW  DY(I) . 
// *** INPUT *** 
// NCD..............SPACE DIMENSION.
// X[NPOINT]........ordinates of the data points.
// DY(I)............ERROR ESTIMATE, 1<=I<=NPOINT. 
// Y(NPOINT,J)......ORRESPONDING DATA ORDINATES at X[i] OF NCD SPACE DIMENSION
//                    (1<=J<=NCD). 
// NPOINT...........NUM OF POINTS. 
// *** OUTPUT *** 
// V(NPOINT,7) 
// QTY(NPOINT,NCD)..Q-TRANSP. * Y
void blgsm1_(int ncd, const double *x, const double *dy, 
	const double *y, int npoint, int iy, double *v,double *qty
){
	double d1, d2, d3;
    double diff, prev;
    int i, j, nm1;
	int n2,n3,n4,n5,n6,n7;

    // Parameter adjustments 
    qty -= npoint+1;
    v -= npoint+1;
    --dy;
    --x;
    y -= iy+1;

    // Function Body 
    nm1=npoint-1;
	n2=npoint*2;
	n3=npoint*3;
	n4=npoint*4;
	n5=npoint*5;
	n6=npoint*6;
	n7=npoint*7;
    v[n4+1] = x[2]-x[1];
    if(npoint<=2)
		return;

    for(i=2; i<=nm1; ++i){
		v[i+n4] = x[i+1]-x[i];//v(i,4)
		v[i+npoint] = dy[i-1]/v[i-1+n4];//v(i,1)
		v[i+n2] = -dy[i]/v[i+n4]-dy[i]/v[i-1+n4];//v(i,2)
		v[i+n3] = dy[i+1]/v[i+n4];//v(i,3)
    }
    v[npoint+npoint] = 0.f;
    for(i=2; i<=nm1; ++i){
		d1 = v[i+npoint];
		d2 = v[i+n2];
		d3 = v[i+n3];
		v[i+n5] = d1*d1+d2*d2+d3*d3;
    }
    if(nm1>=3){
	   for(i=3; i<=nm1; ++i)
			v[i-1+n6] = v[i-1+n2]*v[i+npoint]+v[i-1+n3]*v[i+n2];
    }
    v[nm1+n6] = 0.f;
    if(nm1>=4) {
	    for(i=4; i<=nm1; ++i)
			v[i-2+n7] = v[i-2+n3]*v[i+npoint];
    }
    v[nm1-1+n7] = 0.f;
    v[nm1+n7] = 0.f;
	// CONSTRUCT  Q-TRANSP. * Y  IN  QTY. 
    for(j=1; j<=ncd; ++j){
		prev = (y[j*iy+2]-y[j*iy+1])/v[n4+1];
		for(i=2; i<=nm1; ++i){
		    diff = (y[i+1+j*iy]-y[i+j*iy])/v[i+n4];
		    qty[i+j*npoint] = diff-prev;
		    prev = diff;
		}
    }
}
