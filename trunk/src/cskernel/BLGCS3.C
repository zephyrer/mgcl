/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/Bvltan.h"
#include "cskernel/bkdnp.h"
#include "cskernel/blgcs7.h"

// BLGCS5 IS LOCAL SUB OF BLGCS3, WILL FIND POINT SEQ FOR APPROXIMATE 
// TANGENT, AND GET THE APPROXIMATE TANGENT 
//       ISE....=1 KVALI(JS) IS STARTING POINT, =2 ENDING 
void blgcs5_(int ise, int nvi,const double *taui, 
	const int *kvali, int js,const double *vali, int ivi, int ncd,
	double *work, double *rtan, int *iflag
){
    // Local variables 
    int iend;
    double valw[18], tauw[6];	//valw was [6][3] 
    int i, j;
    int istrt, ii;
    int len;

    // Parameter adjustments 
    --kvali;
    --taui;
    vali -= ivi+1;

    // Function Body 
    len = 1;
    istrt = js;
    iend = js;

    if(ise==1){
		while(len<6 && iend<nvi  && kvali[iend+1] == 0){
			++len;	++iend;
		}
    }else{
		while(len<6 && istrt>1 && kvali[istrt-1] == 0){
			++len;	--istrt;
		}
    }
    for(i=1; i<=len; ++i){
		ii = istrt+i-1;
		tauw[i-1] = taui[ii];
		for(j=1; j<=ncd; ++j)
			valw[i+j*6-7] = vali[ii+j*ivi];
    }
    if(kvali[istrt]==-1){
		for(j=1; j<=ncd; ++j)
		    valw[j*6-6] = vali[istrt-1+j*ivi];
    }
    bkdnp_(&len,tauw,valw,6,ncd,1,bkmax_());
    bvltan_(len,valw,tauw,6,ncd,ise,work,work+70,work+81,rtan,iflag);
}

// INTERNAL SUB OF BLGCS3, STORE DATA IN VAL AND KVAL. 
void blgcs4_(int ncd, int *n, double *tau, 
	double *val, int *kval, double taudt,const double *dat, 
	int ksdt, int iv, int id
){
    int nnew, i;

    // Parameter adjustments 
    --tau;
    --kval;
    val -= iv + 1;
    dat -= id + 1;

    // Function Body 
    nnew = *n + 1;
    tau[nnew] = taudt;
    kval[nnew] = ksdt;
    for(i=1; i<=ncd; ++i)
		val[nnew+i*iv] = dat[i*id+1];
    *n = nnew;
}

// BLGCS3         BLGCS3 FOR BLGCS, RESOLVE KNUCKLE INF KVALI(.). 
// *** INPUT * 
//       IBCI(2)....SPECIFY BOUNDARY COND, WHEN KVALI(1)=0,KVALI(NVI)=0. 
//              =1: 1ST DERIV     =2: 2ND DERIV     =3: NO B.C. 
//              (1ST OR 2ND DERIV. PROVIDED IN TS(.) OR TE(.) ) 
//       TS(NCD),TE(NCD)....NECESSARY ONLY WHEN KVALI(.)=0 AND 
//              IBCI(.)=1 OR 2, AND GIVES 1ST OR 2ND DERIV AS ABOVE. 
//       NCD....SPACE DIMENSION OF THE INPUT POINTS, 1<= NCD <=3. 
//       NVI....NUMBER OF INPUT POINTS, I.E. 
//              ( TAUI(I),KVALI(I),VALI(I,.) )  FOR 1<=I<=NVI 
//       TAUI(NVI),KVALI(NVI),VALI(IVI,NCD)....INPUT POINTS WITH KL INF 
//              ( TAUI(I),KVALI(I),VALI(I,.) ) IS ONE PAIR OF THE INPUT. 
//              TAUI(I) IS A DATA POINT AND VALI(I,.), THE CORRESPND DATA. 
//              KVALI(I) IS THE KNUCKLE INF, AND VALI(.,.) IS THE POINT 
//              (I.E. POSITIONAL DATA OF NCD SPACE DIMENSION) 
//       IVI,IV....ROW DIMENSIONS OF THE VARIABLE VALI AND VAL ,EACH. 
// *** OUTPUT * 
//       IBC(2)....GIVE BOUNDARY CONDITIONS FOR BLG4SQ 
//                 (1)=BEGINNING COND, (2)=ENDING COND. 
//                 =1: 1ST DERIV SPECIFIED 
//                 =2: 2ND DERIV SPECIFIED 
//                 =3: NO BOUNDARY COND. 
//       N,TAU(N),VAL(IV,NCD)....GENERATED DATA SEQUENCE INCLUDING 
//              DERIVATIVE INF. THESE CAN BE STRAIGHT INPUT OF BLG4SQ. 
//                N: NUMBER OF DATA 
//                TAU(I),VAL(I,.):  1<=I<=N ARE DATA POINTS AND 
//                    ASSOCIATED DATA. 
//       IFLAG       INDICATES IF SUCCESSFUL OR NOT, 
//                        = 1 SUCCESSFUL RETURN 
//                        <>1 FAILUE, VAL AND TAU ARE INVALID. 
// ***WORK* 
//   WORK(105),KWORK(N) WORK ARRAY OF EACH LENGTH 
//   STATE NUM.(ABOVE L) ARE AS FOLLOWS 
//      0->0,2->0,1->0 STATE =0           1->1,2->1,0->1 STATE =1 
//      0->2 STATE           =2           1->2,2->2 STATE      =3 
//      1->-1 STATE          =4 
//      X->-1                =5 (X<>1) 
void blgcs3_(const int *ibci,const double *ts,const double *te, 
	int ncd, int nvi,const double *taui,const int *kvali, 
	const double *vali, int ivi, int iv, double *work, 
	int *kwork, int *ibc, int *n, double *tau, double *val, int *iflag
){
    // System generated locals 
    int vali_offset, val_offset;

    // Local variables 
    double told, rtan[3], tnew;
    int i, j, k;
    double r__;
    int kstrt;
    int js, lstate, nm1, nm2, nm3, jsm1;

// **************************** START OF BLGCS3 *************************
    // Parameter adjustments 
    --ibci;
    --te;
    --ts;
    --taui;
    --kvali;
    vali_offset = ivi + 1;
    vali -= vali_offset;
    val_offset = iv + 1;
    val -= val_offset;
    --ibc;
    --tau;

    // Function Body 
    if(nvi<2){
		*iflag = 24;
		return;
    }

// =====START POINT PROCESS===== 
    tnew = taui[1];
    *n = 0;
    blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,tnew,&vali[ivi+1],1,iv,ivi);

    kstrt = kvali[1];
    if (kvali[2] == -1) {
	blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &vali[ivi+2], 1, iv, ivi);
	js = 3;
	ibc[1] = 1;
	lstate = 0;
    } else {
	js = 2;
	if (kstrt == 0) {
	    lstate = 0;
	    if (ibci[1] == 1 || ibci[1] == 2) {
		blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &ts[1], 1, iv, 1);
	    }
//         SET BOUNDARY COND. 
	    ibc[1] = ibci[1];
	} else if (kstrt == 1) {
	    lstate = 1;
	    ibc[1] = 1;
	} else {
	    *iflag = 21;
		*n=js;
		return;
	}
    }
// =====END OF START POINT PROCESS===== 
//  LOOP OVER JS UNTIL JS>NVI 
	while(js<=nvi){
	    told=tnew;
		tnew=taui[js];
		//    GET NEXT K 
	    k = kvali[js];
	    if((lstate==0 || lstate==3 || lstate==4) && k==0){
		// ...CASE OF 0-0,2-0 STATE . <<< LSTATE=0 >>>  NO D.P. MULTIPLICITY. 
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &vali[js+ivi], 0, iv, ivi);
			lstate = 0;
		}else if((lstate==1 || lstate==2) && k==1 || (lstate==1 || lstate==2) && k==2){
		// ...CASE OF 1-1,2-1, 1-2, 2-2 STATE.(STRAIGHT LINE) D.P. MULT=3 
		//                           <<< LSTATE=1 OR 3 >>> 
		//       GET 1ST DERIV. FROM THE STRAIGHT LINE 
			jsm1 = js - 1;
			r__ = taui[js] - taui[jsm1];
			for(i=1; i<=ncd; ++i)
				rtan[i-1] = (vali[js+i*ivi]-vali[jsm1+i*ivi])/r__;
			//       STORE STRAIGHT LINE SPAN 
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, told, rtan, 1, iv, 1);
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, rtan, 1, iv, 1);
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &vali[js+ivi], 1, iv, ivi);
			if(lstate==2){
			    nm2 = *n - 2;
			    nm3 = *n - 3;
				for(i=1; i<=ncd; ++i){
					val[nm2+i*iv] = val[nm3+i*iv];
					val[nm3+i*iv] = rtan[i-1];
				}
			}
			if(k==1)
				lstate = 1;
			else
			    lstate = 3;
	    }else if((lstate==0 || lstate==4) && k == 2){
		// ...CASE OF 0-2 STATE  <<< LSTATE=2 >>>  D.P. MULTIPLICITY = 2. 
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &vali[js+ivi], 1, iv, ivi);
			lstate = 2;
		}else if(lstate==0 && k==1){
		// ...CASE OF 0-1 STATE  <<< LSTATE=1 >>> 
		//     GET APPROXIMATE TANGENT VECTOR ENDING AT THE CURRENT POINT 
			blgcs5_(2,nvi,&taui[1],&kvali[1],js,&vali[vali_offset],ivi,ncd,work,rtan,iflag);
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, rtan, 1, iv, 1);
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tnew, &vali[js+ivi], 1, iv, ivi);
			lstate = 1;
	    }else if(lstate==1 && k==0){
		// ...CASE OF 1-0 STATE  <<< LSTATE=0 >>> 
		//   GET APPROXIMATE TANGENT VECTOR STARTING AT THE PREVIOUS POINT 
			blgcs5_(1,nvi,&taui[1],&kvali[1],js-1,&vali[vali_offset],ivi,ncd,work,rtan,iflag);
			blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,told,rtan,1,iv,1);
			blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,tnew,&vali[js+ivi],0,iv,ivi);
			lstate = 0;
	    }else if(lstate==1 && k==-1){
		// ...CASE OF 1-(-1) STATE <<< LSTATE=4>>> 
			blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,told,&vali[js+ivi],1,iv,ivi);
			lstate = 4;
		}else if((lstate==0 || lstate==3 || lstate==4) && k==-1){
		// ...CASE OF X-(-1) STATE <<< LSTATE=5 >> ( X<>1 ) 
			lstate = 5;
	    }else if(lstate==5 && k==1 || lstate==5 && k==0){
		// ...CASE OF (-1)-1 STATE<< LSTATE=1>>, OR(-1)-0 STATE<< LSTATE=0>> 
			blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,tnew,&vali[js-1+ivi],1,iv,ivi);
			blgcs4_(ncd,n,&tau[1],&val[val_offset],kwork,tnew,&vali[js+ivi],1,iv,ivi);
			lstate=k;
	    }else{
		    *iflag = 22;
			*n=js;
			return;
		}
	    ++js;
	}

// =====END POINT PROCESS=== 
    if(lstate==1 || lstate==0 && kwork[*n-1]==1)
		ibc[2] = 1;
    else if(lstate==0){
		if (ibci[2] == 1 || ibci[2] == 2) {
			blgcs4_(ncd, n, &tau[1], &val[val_offset], kwork, tau[*n], &val[*n + iv], 1, iv, iv);
		    nm1 = *n - 1;
		    for (i = 1; i <= ncd; ++i)
				val[nm1+i*iv] = te[i];
			kwork[nm1-1] = 1;
		}
		ibc[2] = ibci[2];
	}else {
		*iflag = 23;
		*n = js;
		return;
	}
// =====END OF END POINT PROCESS===== 
// SUCCESSFUL RETURN. 
// REMOVE TOO NEAR POINTS 
    blgcs7_(&tau[1], &val[val_offset], kwork, n, ncd, iv);
    if(*n<4 && ibc[1]==3 || *n>=4 && tau[2]==tau[4] && ibc[1]==3){
//   INSERT 2ND DERIV. TO RESCUE ERROR RETURN AT THE START 
		for (i = *n; i >= 2; --i) {
			tau[i + 1] = tau[i];
		    for(j=1; j<=ncd; ++j)
				val[i+1+j*iv] = val[i+j*iv];
		}
		tau[2] = tau[1];
		for(j=1; j<=ncd; ++j)
			val[j*iv+2] = 0.f;
		++(*n);
		ibc[1] = 2;
    }
    if(*n<4 && ibc[2]==3 || *n>=4 && tau[*n-1]==tau[*n-3] && ibc[2]==3){
	//     INSERT 2ND DERIV. TO RESCUE ERROR RETURN AT THE END 
		tau[*n + 1] = tau[*n];
		for(j=1; j<=ncd; ++j)
		    val[*n+1+j*iv] = val[*n+j*iv];
		for(j=1; j<=ncd; ++j)
			val[*n+j*iv] = 0.f;
		++(*n);
		ibc[2] = 2;
    }
    *iflag = 1;
}
