/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/Bvunit.h"
#include "cskernel/bvcsl.h"
#include "cskernel/bvcang.h"

// BLGCS1         BLGCS1 FOR BLGCS 
// BLGCS1 INSERTS CIRCLE DATA INTO GIVEN POINT SEQUENCE VALI(.,.), 
// AND GENERATE NEW KVAL(.,.) 
// *** INPUT *** 
// NCD.....SPACE DIMENSION OF VALI, MUST BE 2 OR 3 WHEN NK4>0. 
// NVI,KVALI(NVI),VALI(IVI,NCD)....INPUT POINT SEQUENCE WITH KNUCKLE 
//         INF.  NVI: NUM OF DATA            KVALI(.): KNUCKLE INF. 
//               VALI(.,NCD): POINTS OF NCD SPACE DIMENSION. 
// NK4,IDK(NK4),RCIR(NK4).....PROVIDE CIRCLE DATA. 
//       RCIR(I) IS RADIOUS OF THE OSCULATING CIRCLE AT THE IMAGINARY 
//       POINT VALI(L,.), I.E. KVALI(L)=1, 
//         WHERE L=IDK(I).       1<=I<=NK4 
//           IDK(.) MUST BE INCREASING, I.E. IDK(I)<IDK(I+1). 
// IVI,IV.....ROW DIMENSION OF VALI AND VAL, EACH. 
// *** OUTPUT *** 
// NV,KVAL(NV),VAL(IV,NCD)......OBTAINED NEW POINT SEQUENCE 
//         INCLUDING CIRCLE DATA. 
//         KVAL AND VAL MAY BE THE SAME AREA AS KVALI AND VALI, EACH. 
// IFLAG....=1: SUCCESSFUL RETURN     <>1:FAILURE 
// *** NOTE * 
//     WHEN KVAL(I)=-1, VAL(I,.) IS NOT POSITIONAL DATA, IS DERIVATIVE. 
//   BLGCS1 GENERATE -1  AS NEW KNUCKLE INF. THESE DATA ARE 
//   RESOLVED BY SUBROUTINE BLGCS3. SEE BLGCS3.
void blgcs1_(int ncd, int nvi, const int *kvali, 
	const double *vali, int nk4, const int *idk,const double *rcir, 
	int ivi, int iv, int *nv, int *kval, double *val, int *iflag
){
    // Initialized data 
    static double themin = .27;

    // Local variables 
    int idki, iaft;
    double tane[3], tans[3];
    int iadd1, iadd2, i, j, n;
    double p[3], r, t, angle;
    double fndiv;
    int nvold, iprev, nipnt;
    double g1[3];
    int i1, i2, i3, nvnew;
    double g2[3], error, p1[3], p2[3], s1, t1, t2, s2;
    double cotang;
    int kesave;
    double center[3];
    int kssave;
    int km1, kp1;
    double ang1, ang2;

    // Parameter adjustments 
    --kval;
    --kvali;
    --rcir;
    --idk;
    vali -= ivi + 1;
    val -= iv + 1;;

    // Function Body 

// ********** START OF THE PROCESS 
// ====== 1. INITIAL SET ====== 
    error = 1.f / bkmax_();
    *iflag = 2;
    for(i=1; i<=nvi; ++i){
		kval[i] = kvali[i];
		for(j=1; j<=ncd; ++j)
			val[i+j*iv] = vali[i+j*ivi];
    }
    kssave=kvali[1];
    if(kssave==7)
		kval[1]=1;
    kesave=kvali[nvi];
    if(kesave==7)
		kval[nvi] = 1;
    nvnew = nvi;
// ======= 2.INSERT CIRCLE DATA ====== 
// LOOP OVER I FROM LAST POINT OF CIRCLE TO 1ST. 
    for(i=nk4; i>=1; --i){
		nvold = nvnew;
		idki = idk[i];
		r = rcir[i];
		iprev = idki - 1;
		iaft = idki + 1;
		km1 = kval[iprev];
		kp1= kval[iaft];
		if (! (km1 >= 1 && km1 <= 3 && (kp1 >= 1 && kp1 <= 3)))
			return;
		//   ... 2.1 COMPUTE OSCULATING CIRCLE. 
		for (j = 1; j <= ncd; ++j) {
		    p[j - 1] = val[idki + j * iv];
			g1[j - 1] = val[iprev + j * iv] - p[j - 1];
		    g2[j - 1] = val[iaft + j * iv] - p[j - 1];
		}
		bvunit_(g1, ncd, tans);
		bvunit_(g2, ncd, tane);
		for(j=1; j<=ncd; ++j)
			tans[j-1] = -tans[j-1];
		//     TANS, TANE ARE UNIT TANGENT VECTORS AT P1,P2 EACH. 
		bvcsl_(r,ncd, p, g1, g2, center, p1, p2, &t1, &t2);
		//     NOW FOLLOWING DATA ARE OBTAINED: 
		//         CENTER(.)=CENTER OF CIRCLE. 
		//         P1(.),P2(.) = START AND END CONTACT POINT. 

	//  ... 2.2 COMPUTE NUMBER OF POINTS TO INSERT. 
		g1[2] = 0.f;
		g2[2] = 0.f;
	//      G1 AND G2 ARE CLEARED FOR THE CASE OF NCD=2. 
		for(j=1; j<=ncd; ++j){
			g1[j-1] = p1[j-1]-center[j-1];
		    g2[j-1] = p2[j-1]-center[j-1];
		}
	//      G1 AND G2 ARE TWO VECTORS OF THE ARC. 
		t = bvcang_(ncd,g1,g2);
		angle = acos(t);
		cotang = t / sin(angle);
		nipnt = (int) (angle / themin);
		if(nipnt<2)
			nipnt=2;

	//  ... 2.3 DETERMINE CONTACT POINTS ADD-KIND. 
	//       IADD=1: ONLY KVAL=2 
	//       IADD=2: KVAL=-1 AND 1 (TANGENT AND CONTACT POINT) 
	//       IADD=3: KVAL=-1 AND 0 
	//     ... 2.3.1 STARTING CONTACT POINT. 
		if (iprev == 1 && kssave == 7)
	    iadd1 = 3;
		else if(fabs(t1-1.f)<=error){
	//         CHECK OF TOO NEAR POINT OF GIVEN AND STARTING CONTACT POINT. 
		    if(km1==2)
				iadd1 = 3;
		    else
				iadd1 = 2;
		}else
		    iadd1 = 1;
	//     ... 2.3.2 ENDING CONTACT POINT. 
		if (iaft == nvi && kssave == 7)
		    iadd2 = 3;
		else if(fabs(t2 - 1.f)<=error){
		//        CHECK OF TOO NEAR POINT OF GIVEN AND ENDING CONTACT POINT. 
		    if (kp1 == 2)
				iadd2 = 3;
			else
				iadd2 = 2;
		}else
		    iadd2 = 1;
		if(iadd1!=1) 
			--iprev;
		if(iadd2!=1)
			++iaft;
	//     ... 2.4 INSERT CIRCLE DATA INTO VAL(.,.) 
		n = nvold-iaft+1;
		nvnew = iprev+2+nipnt+n;
		if(iadd1!=1)
			++nvnew;
		if(iadd2!=1)
			++nvnew;
		i2 = nvnew;
		i3 = nvold;
	//      ... 2.4.1 SHIFT DATA TO GET AREA TO INSERT 
		for(i1=1; i1<=n; ++i1){
		    for(j=1; j<=ncd; ++j)
				val[i2+j*iv] = val[i3+j*iv];
		    kval[i2] = kval[i3];
		    --i2;
			--i3;
		}
	//     ... 2.4.2 INSERT START POINT DATA(TAN VEC AND POINT ATA) 
		n = iprev+1;
		if (iadd1 == 1) {
			kval[n] = 2;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = p1[j-1];
		}else if(iadd1==2){
		    kval[n] = 1;
			for(j=1; j<=ncd; ++j)
				val[n+j*iv] = p1[j-1];
			++n;
		    kval[n] = -1;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = tans[j-1];
		}else{
		    kval[n] = -1;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = tans[j-1];
		    ++n;
			kval[n] = 0;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = p1[j-1];
		}
//      ... 2.4.3 INSERT MID POINTS. 
		fndiv = (double) (nipnt + 1);
		for(i1=1; i1<=nipnt; ++i1){
			t = (double)i1/fndiv;
		    ang1 = t*angle;
		    s1 = cos(ang1)-cotang*sin(ang1);
			ang2 = angle-ang1;
		    s2 = cos(ang2)-cotang*sin(ang2);
		    ++n;
			for(j=1; j<=ncd; ++j)
				val[n+j*iv] = center[j-1]+g1[j-1]*s1+g2[j-1]*s2;
			kval[n] = 0;
		}
	// .... 2.4.4 INSERT END POINT DATA (POINT AND TANGENT) 
		++n;
		if (iadd2 == 1) {
		    kval[n] = 2;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = p2[j-1];
		}else{
		    kval[n] = -1;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = tane[j-1];
		    ++n;
			kval[n] = 1;
		    for(j=1; j<=ncd; ++j)
				val[n+j*iv] = p2[j-1];
		    if(iadd2==3)
				kval[n] = 0;
		}
    }
    *nv = nvnew;
    *iflag = 1;
}
