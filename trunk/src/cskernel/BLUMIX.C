/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/bler.h"
#include "cskernel/Bleval.h"
#include "cskernel/bkdtpg.h"
#include "cskernel/Ble.h"
#include "cskernel/blumi1.h"
#include "cskernel/blumi2.h"
#include "cskernel/bk1fli.h"
#include "cskernel/bkdnp.h"
#include "cskernel/bkmlt.h"
#include "cskernel/blel.h"

// BLUMIX WILL MIX TWO TWO DIMENSIONAL B-REP INTO ONE THREE D B-REP. 
// *** INPUT * 
//  ERROR......ERROR ESTIMATE OF BLIP 
//  KCOD1,K1,N1,T1(N1+K1),RCOEF(IRC1,2).....DESCRIBE THE 1ST B-REP 
//           OF ORDER K1. KCOD1 IS COORDINATE KIND OF THE B-REP, MUST BE 
//           1,2, OR,3 (YZ,ZX,XY EACH) 
//  KCOD2,K2,N2,T2(N2+K2),RCOE2(IRC2,2).....DESCRIBE THE SECOND B-REP OF 
//           ORDER K2. KCOD2 MAY BE 4, INDICATING GIRTH B-REP, I.E. 
//           RCOE2(.,1) IS THE DATA POINTS OF THE 1ST B-REP 
//           AND RCOE2(.,2) IS THE MISSING COORDINATE OF THE 1ST B-REP. 
//  IRC......ROW DIMENSION OF THE VARIABLE RCOEF, I.E. LENGTH OF RCOEF 
//           ALLOWED TO USE. 
//           IRC MUST BE GREATER THAN 11 AND N1. 
// *** OUTPUT * 
//  VECS(3),VECE(3).......START AND END 1ST DERIVATIVES. 
//  NOUT,T(N),RCOEF(IRC,3)......THREE DIMENSIONAL DATA POINTS 
//                           (ORDINATE AND ABSCISSA). 
//  IFLAG........OUTPUT IF SUCCESSFUL RETURN OR NOT. 
//           =0: SUCCESSFUL RETURN 
//           =2: MEMORY EXAUSTED. 
// ***WORK AREA*** 
//  KSEQ(IRC),WORK(NN)....WHERE NN=MAX(IRC*4,4*K*K+3*K).(K=MAX(K1,K2)) 
void blumix_(double error, int kcod1, int k1, 
	int n1,const double *t1,const double *rcoe1, int irc1, 
	int kcod2, int k2, int n2,const double *t2,const double *rcoe2,
	int irc2, int irc, int *kseq, double *work, 
	double *vecs, double *vece, int *nout, double *t, 
	double *rcoef, int *iflag)
{
    static int ksnew[3] = { -2,-1,0 };

	int irc2p1, ircp1, isfm1, ksi;
    int jend, ndel;
    int n1pk1, left;
    double fnew[9];	// was [3][3]
    int mltf, ipos;
    double t1old, t2old, f[3], g[3];
    int n, i, j, l, kcod11, kcod12;
    int kcodx;
    int i1, j1, j2, jstrt;
    double t1para, t2para;
    int kk;
    double taunew[3];
    int k1m2, k2m1, n1p1, n2p1, isf, nin;
    double tin, tan1[3], tan2[3], tin1;

    // Parameter adjustments 
    --t1;
    irc2p1 = irc2 + 1;
    rcoe2 -= irc2p1;
    --t2;
    ircp1 = irc + 1;
    rcoef -= ircp1;
    --t;
    work -= ircp1;
    --kseq;
    --vece;

    // Function Body 

    kcod11 = kcod1%3+1;
    kcod12 = kcod11%3+1;
    if(kcod2==4){
		kcodx = 2;
    }else if(kcod11==kcod2){
		kcodx = 2;
    }else{
		kcodx = 1;
    }
    n1p1 = n1+1;
    n1pk1 = n1+k1;
    n2p1 = n2+1;
    k1m2 = k1-2;
    k2m1 = k2-1;

// === GET DATA OF B-REP1,2 , DATA KIND KSEQ(.), AND DATA POINT T(.) ===
    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1[k1],1,tan1);
    ble_(k2,n2,&t2[1],&rcoe2[irc2p1], irc2, 2, t2[k2], 1, tan2);
    blumi4_(kcod1, tan1, kcod2, tan2, vecs, 1);
    t[1] = t1[k1];
    kseq[1] = 0;
    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1[k1],0,f);
    rcoef[kcod11*irc+1] = f[0];
    rcoef[kcod12*irc+1] = f[1];
    rcoef[kcod1*irc+1] =
		bler_(k2,n2,&t2[1],&rcoe2[kcodx*irc2+1],t2[k2],0);
    n = 2;
    i1 = k1 + 1;
    t2old = t2[k2];

    while(i1<=n1p1){
		bkmlt_(n1pk1, &t1[1], i1, k1m2, &isf, &mltf);
		if(isf==-1 || isf>n1)
			isf = n1p1;

		//   .....(2.1) INSERT DATA OF CONTINUOUS POINTS. 
		isfm1 = isf-1;
		for(i=i1; i<=isfm1; ++i){
			kseq[n] = 0;
		    t[n] = t1[i];
		    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1[i],0,f);
			rcoef[n+kcod11*irc] = f[0];
		    rcoef[n+kcod12*irc] = f[1];
		    blumi2_(kcod2,k2,n2,&t2[1],&rcoe2[irc2p1],irc2,kcod1, 
			    t[n],f,t2old,error,&work[ircp1],&t2para);
		    rcoef[n+kcod1*irc] =
				bler_(k2,n2,&t2[1],&rcoe2[kcodx*irc2+1],t2para,0);
			t2old = t2para;
		    ++(n);
		}

		//   .....(2.2) INSERT DATA OF DISCONTINUOUS POINTS. 
		i1 = isf;
		if(mltf==k1m2 || i1==n1p1){
			kseq[n] = -1;
		    kseq[n+1] = 0;
		    t[n] = t1[i1];
			t[n+1] = t1[i1];
		}else{
		    kseq[n] = -2;
			kseq[n+1] = -1;
		    kseq[n+2] = 0;
		    t[n] = t1[i1];
			t[n+1] = t1[i1];
		    t[n+2] = t1[i1];
		}

		ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1[i1],0,f);
		rcoef[n+1+kcod11*irc] = f[0];
		rcoef[n+1+kcod12*irc] = f[1];
		blumi2_(kcod2,k2,n2,&t2[1],&rcoe2[irc2p1],irc2,kcod1,t[n],
			f,t2old,error,&work[ircp1],&t2para);
		left=bk1fli_(n2p1,&t2[1],t2para);
		if(t2para-t2[left]<=error){
		    t2para = t2[left];
		}else if(t2[left+1]-t2para<=error)
		    t2para = t2[left+1];

		//      STORE LEFT CONTINUOUS 1ST DERIV. 
		bleval_(k1,n1,&t1[1],rcoe1,irc1,2,1,t[n],1,2,1,tan1);
		bleval_(k2,n2,&t2[1],&rcoe2[irc2p1],irc2,2,1,t2para,1,2,1,tan2);
		blumi4_(kcod1,tan1,kcod2,tan2,&rcoef[n+irc],irc);
		++n;
		rcoef[n+kcod1*irc] =bler_(k2,n2,&t2[1],&rcoe2[kcodx*irc2+1],t2para, 0);
		++n;
		if(mltf!=k1m2 && i1<n1p1) {
		//         STORE RIGHT CONTINUOUS 1ST DERIV. 
		    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t[n],1, tan1);
		    ble_(k2,n2,&t2[1],&rcoe2[irc2p1],irc2,2,t2para,1,tan2);
			blumi4_(kcod1, tan1, kcod2, tan2, &rcoef[n + irc], irc);
		    ++n;
		}
		t2old = t2para;
		i1 += mltf;
    }
    --n;
    for(j=1; j<=3; ++j){
		vece[j] = rcoef[n-1+j*irc];
		rcoef[n-1+j*irc] = rcoef[n+j*irc];
    }
    --n;
    kseq[n] = 0;

// NOW 3-D DATA OF BREP1 KNOT ARE OBTAINED. 
// === (3). INSERT ADDITIONAL DATA OF B-REP2 (IF NECESSARY). === 
    i = 1;
    jstrt = k2;
    t1para = t[i];
    t2old = t2[k2];
    isf = 0;
    while(isf!=-1 && i<n){
		bkmlt_(n2, &t2[1], jstrt, k2m1, &isf, &mltf);
		jend = isf;
		if(isf==-1)
		    jend = n2;
		++i;
		j1 = jstrt;
		t1old = t1para;
		//  ....(3.1) BREP2 DATA INSERT AT CONTINUOUS POINTS 
		while(j1+1<jend && i<n){
		    ipos = i - 1;
			if (kseq[i] < 0) {
			++i;
		    }
			f[0] = rcoef[i + kcod11 * irc];
		    f[1] = rcoef[i + kcod12 * irc];
		    if (kseq[i - 1] == -2) {
			++i;
		    }
		    blumi2_(kcod2,k2,n2,&t2[1],&rcoe2[irc2p1],irc2,kcod1, 
			    t[i],f,t2old,error,&work[ircp1],&t2para);
		    j2=bk1fli_(n2p1,&t2[1],t2para);
		    if(j2>jend){
				j2 = jend;
				t2para = t2[j2];
		    }
		    nin = j2-j1-1;
			//      .....When I jumped too many knots of brep2,insert data of B-rep2. 
		    for(l=1; l<=nin; ++l){
				tin = t2[j1+l-1];
				if(l==1)
				    tin = t2old;
				tin = tin+t2[j1+l]+t2[j1+l+1];
				tin1 = t2[j1+l+2];
				if(l==nin)
					tin1 = t2para;
				tin = (tin+tin1)/4.f;
				ble_(k2,n2,&t2[1],&rcoe2[irc2p1],irc2,2,tin,0,f);
				blumi2_(kcod1,k1,n1,&t1[1],rcoe1,irc1, 
					kcod2,tin,f,t1old,error,&work[ircp1],&t1para);
				f[kcod1-1] = f[kcodx-1];
				ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1para,0,g);
				f[kcod11-1] = g[0];
				f[kcod12-1] = g[1];
				//           Now data to insert in F(1-3) 
				taunew[0] = t1para;
				blumi1_(&n,&kseq[1],&rcoef[ircp1],irc,&t[1],ipos, 
					1,&ksnew[2],f,taunew,iflag);
				if(*iflag!=0)
				    return;
				t1old = t1para;
				++ipos;
				++i;
		    }
		    t1old = t[i];
		    t2old = t2para;
			++i;
		    j1 = j2;
		}
	//  ....  (3.2) Insert data of discontinuous Point of Brep-2. 
		if(isf!=-1){
			t2old = t2[isf];
		    jstrt = jend+mltf;
			//		    Get associated point of Brep-2 
			ble_(k2,n2,&t2[1],&rcoe2[irc2p1],irc2,2,t2[isf],0,f);
		    blumi2_(kcod1,k1,n1,&t1[1],rcoe1,irc1,
				kcod2,t2[isf],f,t1old,error,&work[ircp1],&t1para);
			//      Check if corresponding point in T. 
		    i=bk1fli_(n, &t[1], t1para);
			if(t1para-t[i]<=error){
				ipos = i;
				if(kseq[ipos-1]<0)
					--ipos;
				if(kseq[ipos-1]<0)
					--ipos;
				t1para = t[i];
				if(kseq[ipos]== -2)
					continue;
				ndel = 1 - kseq[ipos];
		    }else if(t[i+1]-t1para<=error){
				++i;
				ipos = i;
				t1para = t[i];
				kk = kseq[i];
				i -= kk;
				if(kk==-2)
					continue;
				ndel = 1 - kk;
		    }else{
				ndel = 0;
				ipos = i;
		    }
			//    ...Delete data of Brep1 from RCOEF 
		    if(ndel>0){
				blumi3_(&n,&kseq[1],&rcoef[ircp1],irc,&t[1],ipos,ndel);
				--ipos;
		    }
			//   ... Get three D data from Brep-1 an 2. 
			tan1[0] = blel_(k1,n1,&t1[1],rcoe1,t1para,1);
		    tan1[1] = blel_(k1,n1,&t1[1],rcoe1+irc1+1,t1para,1);
		    tan2[0] = blel_(k2,n2,&t2[1],&rcoe2[irc2+1],t2[isf],1);
			tan2[1] = blel_(k2,n2,&t2[1],&rcoe2[(irc2<<1)+1],t2[isf],1);
		    blumi4_(kcod1, tan1, kcod2, tan2, fnew, 1);
		    taunew[0] = t1para;
			fnew[kcod1+2] = f[kcodx-1];
		    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1para,0,g);
		    fnew[kcod11+2] = g[0];
			fnew[kcod12+2] = g[1];
		    taunew[1] = t1para;
		    ble_(k1,n1,&t1[1],rcoe1,irc1,2,t1para,1,tan1);
			ble_(k2,n2,&t2[1],&rcoe2[irc2p1],irc2,2,t2[isf],1,tan2);
		    blumi4_(kcod1,tan1,kcod2,tan2,&fnew[6],1);
		    taunew[2] = t1para;
			blumi1_(&n, &kseq[1], &rcoef[ircp1], irc, &t[1], ipos, 3,
				ksnew, fnew, taunew, iflag);
		    if(*iflag != 0)
				return;
			i = ipos + 3;
		}
    }
// ===== (4) Now all necessary data obtained in RCOEF(.,.), get 3d B-rep. 
    l = 0;
    i = 1;
    while(i <= n) {
	mltf = 1 - kseq[i];
	ipos = i;
	if (mltf > 1) {
	    ++ipos;
	}
	++l;
	kseq[l] = mltf;
	for (j = 1; j <= 3; ++j) {
	    work[l+j*irc] = rcoef[ipos+j*irc];
	}
	i += mltf;
    }
    bkdtpg_(&work[ircp1], l, 3, irc, &work[(irc<<2)+1]);
    n = 0;
    for(i=1; i<=l; ++i){
		ksi = kseq[i];
		for(j=1; j<=ksi; ++j){
			++(n);
		    t[n] = work[i+(irc<<2)];
		}
    }
    bkdnp_(&n,&t[1],&rcoef[ircp1],irc,3,2,bkmax_());
	*nout=n;
}
