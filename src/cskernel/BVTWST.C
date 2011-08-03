/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bvabs.h"
#include "cskernel/bz3sol.h"
#include "cskernel/bvvrt.h"
#include "cskernel/bvs3pd.h"
#include "cskernel/Bvltan.h"
#include "cskernel/bvinpd.h"

//  BVTWST COMPUTES TWIST DATA FROM SURFACE DATA SURF(.,.,.), AND 
//  STORES BACK IN SURF. 
// ******  I N P U T  ****** 
//  IBC(4)....SPECIFIES BOUNDARY COND ALONG THE 4 BOUNDARIES, 
//        IBC(I)=1: TANGENT VECTOR ALONG THE PERIMETER I. 
//        IBC(I)<>1: NO B.C. FOR THE PERIMETER I. 
//        ID I IS PERIMETER NUM. OF GIVEN PATCH. 
//        WHEN IBC(I) <> 1, FOLLOWING E(.,I), DEDU(.,I), DEDV(.,I) 
//        ARE INVALID. 
//  E(3,4).....ARE NORMAL VECTORS OF THE SURFACE AT THE FOUR CORNER OF 
//        THE SURFACE. 
//  DEDU(3,4), DEDV(3,4)....ARE DERIVATIVES OF E(,.) ALONG U(DEDU) 
//         AND ALONG V (DEDV). 
//         I=1: U-MIN AND V-MIN, I=2: U-MAX AND V-MIN, 
//          =3: U-MAX AND V-MAX,  =3: U-MIN AND V-MAX 
//  LUD,UTAU(LUD),LVD,VTAU(LVD),SURF(ISR1,ISR2,3)....... 
//        TANGENT DATA INPUT AS FOLLOWS, ACCORDING TO IBC(I); 
//     SURF(I,2,.)     1<=I<=LUD ALONG V=MIN PARAM LINE(TAN VEC OF V DIR). 
//     SURF(LUD-1,J,.) 1<=J<=LVD ALONG U=MAX PARAM LINE(TAN VEC OF U DIR). 
//     SURF(I,LVD-1,.) 1<=I<=LUD ALONG V=MAX PARAM LINE(TAN VEC OF V DIR). 
//     SURF(2,J,.)     1<=J<=LVD ALONG U=MIN PARAM LINE(TAN VEC OF U DIR). 
//  ISR1,ISR2....ARE ARRAY LENGTH OF SURF (ISR1,ISR2) 
// ***** W O R K ***** 
//  WK1(70),WK2(MM)  WORKARRAYS OF EACH LENGTH, WHERE MM=MAX(LUD,LVD) 
// *****  O U T P U T  ***** 
// SURF...CONTAINS TWIST DATA AT THEIR FOUR CORNER; 
//     SURF(2,2,.)          AT U=MIN,V=MIN. 
//     SURF(LUD-1,2,.)      AT U=MAX,V=MIN. 
//     SURF(LUD-1,LVD-1,.)  AT U=MAX,V=MAX. 
//     SURF(2,LVD-1,.)      AT U=MIN,V=MAX. 
// IFLAG....=1 SUCCESSFUL RETURN, 
//          =85 SET AS TWIST DATA = ALL 0.0 
//          <>1 SOME ERROR DETECTED. 
// ***** N O T E ***** 
// ORDER OF ALL B-REP ASSUMED TO BE 4. 
void bvtwst_(const int *ibc,const  double *e,const double *dedu, 
	const double *dedv, int lud,const double *utau, int lvd,const double *vtau,
	double *surf, int isr1, int isr2, 
	double *wk1, double *wk2, int *iflag)
{
    // Local variables 
    double df1e2[4], dsdu[4], wk3[11], wk4[8];
    double dsdv[4], vmin;
    double df1[12], df2[12];	// df1 was [3][4] 	// df2 was [3][4] ;
    double r, a[9], d[3];	//a was [3][3] ;
    double twist[12];	// was [3][4]
    int iflg1,ludm1,lvdm1,i,j,m[4],i1,i2,j1,j2,ni,nj,ntw,idf1[4],idf2[4],isr12;

    // Parameter adjustments 
    --ibc;
    e -= 4;
    dedu -= 4;
    dedv -= 4;
    --wk2;
    --utau;
    --vtau;
    surf -= isr1*(isr2+1)+1;;

    // Function Body 
    vmin = bzmzro_();
    if(lud<=3 || lvd<=3){
		*iflag = 84;
		return;
    }
    *iflag = 1;
    isr12 = isr1*isr2;
    ludm1 = lud-1;
    lvdm1 = lvd-1;

    ntw = 0;
    if(ibc[4]==1 && ibc[1]==1){
		ntw = 1;
		m[0] = 1;
    }
    for(i=2; i<=4; ++i){
		if(ibc[i-1]==1 && ibc[i]==1){
		    ++ntw;
			m[ntw-1] = i;
		}
    }
    if(ntw==0)
		return;

// ====== 1. OBTAIN D(ABS(DS/DU))/DV, etc 
    i1 = 3;
    if(ibc[4]!=1)
		i1 = 2;
    i2 = lud-2;
    if(ibc[2]!=1)
		i2 = ludm1;
    j1 = 3;
    if(ibc[1]!=1)
		j1 = 2;
    j2 = lvd - 2;
    if(ibc[3]!=1)
		j2 = lvdm1;

//     1.1 V=MIN LINE 
    if(ibc[1]==1){
		wk2[1] = bvabs_(3,isr12,&surf[(isr2+2)*isr1+1]);
		ni = 1;
		for(i=i1; i<=i2; ++i){
			++ni;
		    wk2[ni] = bvabs_(3,isr12,&surf[i+(isr2+2)*isr1]);
		}
		++ni;
		wk2[ni] = bvabs_(3,isr12,&surf[lud+(isr2+2)*isr1]);
		bvltan_(ni,&wk2[1],&utau[i1-1],1,1,1,wk1,wk3,wk4,dsdv,iflag);
		bvltan_(ni,&wk2[1],&utau[i1-1],1,1,2,wk1,wk3,wk4,&dsdv[1],iflag);    
	}
//     1.2 V=MAX LINE 
    if(ibc[3]==1){
		wk2[1] = bvabs_(3,isr12,&surf[(lvdm1+isr2)*isr1+1]);
		ni = 1;
		for(i=i1; i<=i2; ++i){
			++ni;
		    wk2[ni] = bvabs_(3,isr12,&surf[i+(lvdm1+isr2)*isr1]);
		}
		++ni;
		wk2[ni] = bvabs_(3, isr12, &surf[lud+(lvdm1+isr2)*isr1]);
		bvltan_(ni,&wk2[1],&utau[i1-1],1,1,1,wk1,wk3,wk4,&dsdv[3],iflag);
		bvltan_(ni,&wk2[1],&utau[i1-1],1,1,2,wk1,wk3,wk4,&dsdv[2],iflag);
    }
//     1.3 U=MIN LINE 
    if(ibc[4]==1){
		wk2[1] = bvabs_(3,isr12,&surf[(isr2+1)*isr1+2]);
		nj = 1;
		for(j=j1; j<=j2; ++j){
		    ++nj;
			wk2[nj] = bvabs_(3,isr12,&surf[(j+isr2)*isr1+2]);
		}
		++nj;
		wk2[nj] = bvabs_(3,isr12,&surf[(lvd+isr2)*isr1+2]);
		bvltan_(nj,&wk2[1],&vtau[j1-1],1,1,1,wk1,wk3,wk4,dsdu,iflag);
		bvltan_(nj,&wk2[1],&vtau[j1-1],1,1,2,wk1,wk3,wk4,&dsdu[3],iflag);
    }
//     1.4 U=MAX LINE 
    if(ibc[2]==1){
		wk2[1]=bvabs_(3,isr12,&surf[ludm1+(isr2+1)*isr1]);
		nj = 1;
		for(j=j1; j<=j2; ++j){
		    ++nj;
			wk2[nj] = bvabs_(3,isr12,&surf[ludm1+(j+isr2)*isr1]);
		}
		++nj;
		wk2[nj] = bvabs_(3,isr12,&surf[ludm1+(lvd+isr2)*isr1]);
		bvltan_(nj,&wk2[1],&vtau[j1-1],1,1,1,wk1,wk3,wk4,&dsdu[1],iflag);
		bvltan_(nj,&wk2[1],&vtau[j1-1],1,1,2,wk1,wk3,wk4,&dsdu[2],iflag);
    }

// ===== 2. OBTAIN (DF1*DE2+DF2*DE1)/2 ( MAY BE INVALID). 
    for(j=1; j<=3; ++j){
		df1[j-1] = surf[(j*isr2+1)*isr1+2];
		df1[j+2] = surf[ludm1+(j*isr2+1)*isr1];
		df1[j+5] = surf[ludm1+(lvd+j*isr2)*isr1];
		df1[j+8] = surf[(lvd+j*isr2)*isr1+2];

		df2[j-1] = surf[(j*isr2+2)*isr1+1];
		df2[j+2] = surf[lud+(j*isr2+2)*isr1];
		df2[j+5] = surf[lud+(lvdm1+j*isr2)*isr1];
		df2[j+8] = surf[(lvdm1+j*isr2)*isr1+1];
    }
    for(i=1; i<=ntw; ++i){
		df1e2[m[i-1]-1] = bvinpd_(3,&df2[m[i-1]*3-3],&dedu[m[i-1]*3+1]);
		df1e2[m[i-1]-1] = (df1e2[m[i-1]-1]+bvinpd_(3,&df1[m[i-1]*3-3],&dedv[m[i-1]*3+1]))* -.45f;
    }

// ===== 3. OBTAIN UNIT OF DF1(.,.) (U-DIR) AND DF2(.,.) (V-DIR) 
// CONVERT DF1,DF2 TO UNIT VECTOR. 
    for(i=1; i<=ntw; ++i){
		r = bvabs_(3,1,&df1[m[i-1]*3-3]);
		if(r<=vmin){
			idf1[i-1] = -1;
		}else{
		    for(j=1; j<=3; ++j)
			df1[j+m[i-1]*3-4] /= r;
		    idf1[i-1] = 0;
		}
		r = bvabs_(3,1,&df2[m[i-1]*3-3]);
		if(r<=vmin){
			idf2[i-1] = -1;
		}else{
		    for(j=1; j<=3; ++j)
				df2[j+m[i-1]*3-4]/= r;
			idf2[i-1] = 0;
		}
    }
// CHECK IF DF WAS ZERO. 
    for(i=1; i<=ntw; ++i){
	if(idf1[i-1]==-1){
	    if(m[i-1]<= 2){
			for(j=1; j<=3; ++j)
				d[j-1]=surf[lud+(lvd+j*isr2)*isr1]-surf[(lvd+j*isr2)*isr1+1];
	    }else{
			for(j=1; j<=3; ++j)
				d[j-1]=surf[lud+(j*isr2+1)*isr1]-surf[(j*isr2+1)*isr1+1];
	    }
	    bvvrt_(&e[m[i-1]*3+1],d,&df1[m[i-1]*3-3]);
	}
    }
    for(i=1; i<=ntw; ++i){
	if(idf2[i-1]==-1){
	    if(m[i-1]==1 || m[i-1]==4){
			for(j=1; j<=3; ++j)
				d[j-1] = surf[(lvd+j*isr2)*isr1+1]-surf[(j*isr2+1)*isr1+1];
	    }else{
			for(j=1; j<=3; ++j)
				d[j-1] = surf[lud+(lvd+j*isr2)*isr1]-surf[lud+(j*isr2+1)*isr1];
	    }
	    bvvrt_(&e[m[i-1]*3+1],d,&df2[m[i-1]*3-3]);
	}
    }

// ===== 4. OBTAIN TWIST VECTORS BY SOLVING LINEAR EQUATIONS. 
    *iflag = 1;
    for(i=1; i<=ntw; ++i){
		d[0] = df1e2[m[i-1]-1];
		d[1] = dsdu[m[i-1]-1];
		d[2] = dsdv[m[i-1]-1];
		for(j=1; j<=3; ++j){
			a[j*3-3] = e[j+m[i-1]*3];
		    a[j*3-2] = df1[j+m[i-1]*3-4];
		    a[j*3-1] = df2[j+m[i-1]*3-4];
		}
		r = bvs3pd_(a,&a[3],&a[6]);
		if(fabs(r)>vmin){
		    bz3sol_(a,d,&twist[m[i-1]*3-3],&iflg1);
		}else{
		    for(j=1; j<=3; ++j)
				twist[j+m[i-1]*3-4] = 0.f;
		    *iflag = 85;
		}
    }
// STORE TWIST DATA. 
    for(j=1; j<=3; ++j){
		if (ibc[4]==1 && ibc[1]==1)
		    surf[(j*isr2+2)*isr1+2] = twist[j-1];
		if (ibc[1] == 1 && ibc[2] == 1) 
			surf[ludm1+(j*isr2+2)*isr1] = twist[j+2];
		if (ibc[2]==1 && ibc[3]==1)
		    surf[ludm1+(lvdm1+j*isr2)*isr1] = twist[j+5];
		if(ibc[3]==1 && ibc[4]==1)
			surf[(lvdm1+j*isr2)*isr1+2] = twist[j+8];
    }
}
