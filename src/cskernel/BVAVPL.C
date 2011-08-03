/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bvvpd.h"
#include "cskernel/bvinpd.h"
#include "cskernel/bvdpsl.h"
#include "cskernel/bvabs.h"

// BVAVPL OBTAINS NEAREST PLANE TO THE GIVEN PLYGON POLY(.,.,.) 
// *** INPUT  * 
//   ERRORI....ERROR ALLOWED TO REGARD A SAME POINT OR SAME S LINE. 
//   LUD,LVD,POLY(IPU,IPV,NCDI),NCDI,IPU,IPV 
//                                .....GIVE POLYGON ARRAY. 
//   (I,J)-TH POINT OF SPACE DIMENSION NCDI IS POLY(I,J,K) (1<=K<=NCDI) 
//   NCDI<=3 IS ASSUMED. 
// *** OUTPUT * 
//   KPLANE, G(4)....TELLS KIND OF PLANE, I.E.: 
//      =0:LUD<=0 && LVD<=0 
//      =1:POLY ARE ON ONE POINT :G(J)  (1<=J<=NCDI) 
//      =2:POLY ARE ON A STRAIGHT LINE; 
//          G(1 to NCDI): DIRECTIONAL COSINE 
// 		 Point on the line is obtained from CNTR(.). 
//      =3: A PLANE G(I) (1<=I<=NCDI+1): G(1)*X+G(2)*Y+G(3)*Z=G(4) 
//   CNTR(NCDI).....CENTER POINT OF POLY(.,..,.) 
//   DEVMAX......IS MAXIMUM DEVIATION OF POLY FROM THE PLANE G(.) 
void bvavpl_(double errori, int lud, int lvd, 
	const double *poly, int ncdi, int ipu, int ipv,
	int *kplane, double *g, double *cntr, double *devmax)
{
    double d;
    double flud, flvd, vmin, ersq;
    int ipuv;
    int i, j, k;
    int isave, jsave;
    double c1, error;
    int i1, i2;
    int ncd;
    double dev;

    double a[6], e[3];// a was [3][2] 
    double p[3], r;
    double plane[3], esave[3],gs[6];

    double el, vl[2];

    // Parameter adjustments 
    --cntr;
    poly -= ipu*(ipv+1)+1;;
    --g;

    // Function Body 
    vmin = bzmzro_();
    *kplane = 0;
    *devmax = 0.f;
    if(lud<1 || lvd<1)
		return;

    error = errori;
    if(error<vmin)
		error = vmin;
    ersq = error*error;
    ncd = ncdi;
    if(ncd>3)
		ncd=3;

// ===== 1. INITIAL SET(GET AVERAGE POINT & CLEAR) ======== 
    flud = (double) (lud);
    flvd = (double) (lvd);
    for(k=1; k<=ncd; ++k){
		cntr[k] = 0.f;
		for(j=1; j<=lvd; ++j){
		    c1 = 0.f;
			for(i=1; i<=lud; ++i)
				c1 += poly[i+(j+k*ipv)*ipu];
		    c1 /= flud;
		    cntr[k] += c1;
		}
		cntr[k] /= flvd;
		//      CNTR(.) IS CENTER (AVERAGE) POINT OF INPUT POLYGON ARRAY. 
		gs[k-1] = cntr[k];
    }
// ===== 2. TEST IF ONE POINT OR NOT. ===== 
    for(j=1; j<=lvd; ++j){
		for(i=1; i<=lud; ++i){
		    dev = 0.f;
		    for(k=1; k<=ncd; ++k){
			// Computing 2nd power 
				d=cntr[k]-poly[i+(j+k*ipv)*ipu];
				dev+=d*d;
			}
		    if(dev>*devmax){
				isave = i;
				jsave = j;
				*devmax = dev;
				if(*devmax>ersq)
					goto L50;
		    }
		}
    }
    *devmax = sqrt(*devmax);
    *kplane = 1;
    return;

L50:
// ===== 3. TEST IF STRAIGHT LINE OR NOT. 
    for(k=1; k<=ncd; ++k)
		e[k-1] = poly[isave+(jsave+k*ipv)*ipu]-cntr[k];
    r = bvabs_(ncd,1,e);
    for(k=1; k<=ncd; ++k)
		gs[k+ncd-1] = e[k-1]/r;

    *devmax = 0.f;
    for(i=1; i<=lud; ++i){
		for(j=1; j<=lvd; ++j){
		    for(k=1; k<=ncd; ++k)
				p[k-1] = poly[i+(j+k*ipv)*ipu];
			dev = fabs(bvdpsl_(ncd,gs,p));
		    if(dev>error)
				goto L150;
			if(*devmax<dev)
				*devmax=dev;
		}
    }
    *kplane=2;
    for(k=1; k<=ncd; ++k)
		g[k] = gs[k+ncd-1];
    return;
// ===== 4. SET INITIAL VECTOR 
L150:
    *kplane = 3;
	if(ncd==2){
		g[1]=0.; g[2]=0.; g[3]=1.; g[4]=0.;
		*devmax=0.;
		return;
	}

    for(k=1; k<=3; ++k){
		plane[k-1] = 0.f;
		esave[k-1] = 0.f;
		a[k-1] = 0.f;
		a[k+2] = 0.f;
    }
    for(k=1; k<=ncd; ++k)
		a[k-1] = poly[(k*ipv+1)*ipu+1]-cntr[k];
    vl[0] = bvabs_(ncd,1,a);
    if(vl[0]>vmin)
		for(k=1;k<=ncd;++k)
			a[k-1]/=vl[0];
// ===== 5. LOOP FROM SECOND VECTOR TO N-TH VECTOR. ===== 
    i1 = 0;
    for(j=1; j<=lvd; ++j){
	for(i=1; i<=lud; ++i){
	    i1 = i1 % 2 + 1;
	    i2 = i1 % 2 + 1;
	    for(k=1; k<=ncd; ++k)
			a[k+i2*3-4] = poly[i+(j+k*ipv)*ipu]-cntr[k];
			vl[i2-1] = bvabs_(ncd,1,&a[i2*3-3]);
		if(vl[i2-1]>vmin){
			for(k=1; k<=ncd; ++k)
				a[k+i2*3-4] /= vl[i2-1];
			if(vl[i1-1]>vmin){
				bvvpd_(&a[i1*3-3],&a[i2*3-3],e);
				el = bvabs_(3,1,e);
				if(el>vmin){
					for (k = 1; k <= 3; ++k) e[k - 1] /= el;
					if(e[0]*esave[0]+e[1]*esave[1]+e[2]*esave[2]<0.f)
						for(k=1; k<=3; ++k)
							e[k-1] = -e[k-1];
					for(k=1; k<=3; ++k){
					    esave[k-1] = e[k-1];
					    plane[k-1] += esave[k-1];
					}
				}
			}
		}
	}
    }

// ===== 4. DEFINE KPLANE AND GET MAXIMUM DEVIATION, DEVMAX ===== 
    el = bvabs_(3,1,plane);
    for(k=1; k<=3; ++k)
		g[k]=plane[k-1]/el;
    g[4] = bvinpd_(ncd,&g[1],&cntr[1]);
//      GET MAXIMUM DEVIATION FORM THE PLANE 
    *devmax = 0.f;
    ipuv = ipu*ipv;
    for(j=1; j<=lvd; ++j){
	for(i=1; i<=lud; ++i){
	    for(k=1; k<=ncd; ++k)
			p[k-1] = poly[i+(j+k*ipv)*ipu];
	    dev = bvinpd_(ncd,&g[1],p)-g[4];
	    if(*devmax<fabs(dev))
			*devmax=fabs(dev);
	}
    }
}
