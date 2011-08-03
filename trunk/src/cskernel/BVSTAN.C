/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/Ble.h"
#include "cskernel/Bleval.h"
#include "cskernel/bvvrt.h"
#include "cskernel/bvivec.h"
#include "cskernel/bvabs.h"
#include "cskernel/blgint.h"
#include "cskernel/bkdtkt.h"

//         SUBROUTINE TO GET 1ST DERIVATIVE FROM THE TANGENT PLANE. 
// ***INPUT* 
//  SET(2).....PROVIDES PARAM VALUE OF THE INTERESTING PERIMETER 
//             SET(1)=STARTING PARAM,     (2)=ENDING 
//  KTP,NTP,TTP(.),RTP(IRTP,3).... PROVIDES B-REP OF TANGENT PLANE (ORDER KTP) 
//                 (UNIT NORMAL VECTOR OF TANGENT PLANE) 
//  TPTAU...PARAMETER VALUE OF THE BOUNDARY AT WHICH P SHOULD BE TANGNT 
//  N,TAU(N),P(IP1,IP2,3).....POSITIONAL DATA SEQ. OF THE LINE 
//           DATA PROVIDED IN ROWWISE, I.E. P(1,I,.) 1<=I<=N 
//           (THE 1ST OR LAST 6 DATA ARE USED FOR APPROXIMATE TANGN) 
//  IPSE......INDICATES WHICH BOUNDARY POINT'S TANGENT BE OBTAINED 
//                STARTING POINT(=1)  OR ENDDING POINT (=2) OF P(.,.) 
//  ITANSE(2)..INDICATE WHETHER TANSE(.,.) ARE PROVIDED OR NOT. 
//            ITANSE(I)=1: TANSE(.,I) PROVIDED,  =2: NOT PROVIDED 
//  TANSE(3,2)....TANSE(.,I) IS VALID ONLY WHEN ITANSE(I)=1, AND 
//            GIVE TANGENT VECTOR OF STARTING OR END POINT OF THE 
//            TANGENT PLANE B-REP. I=1: STARTING,   =2: ENDING 
//            TANSE(.,I) RESTRICT OBTAINING TANGNT. 
//  IRTP,IP1,IP2......ARE ROW DIMENSION OF RTP AND P, EACH. 
// ***OUTPUT* 
//    TANGN(3)....1ST DERIVATIVE AT THE POINT TPTAU WHOSE DATA POINTT 
//            SEQ. ARE TAU(.). 
// ***WORK * 
//    WORK(105)....WORK ARRAY. 
void bvstan_(const double set[2], int ktp, int ntp, 
	const double *ttp,const double *rtp, double tptau, int n, 
	const double *tau,const double *p, int ipse, int *itanse, 
	const double *tanse, int irtp, int ip1, int ip2, 
	double *work, double *tangn){
    // Initialized data 

    static const int ip = 6;
    static const int k = 4;
    static const int ncd = 3;

    // System generated locals 
    int rtp_offset;

    // Local variables 
    double rabs;
    int i, j, iflag;
    double s, r, absse=1.;
    double tlenh, vlmin, point[18];	// was [6][3]
    int istrt;
    double t1, t2;
    int i12, kp;
    int np2;
    double dva[3], dvb[3], dve[3], tpc;

    // Parameter adjustments 
    --ttp;
    --tau;
    --itanse;
    tanse -= 4;
    rtp_offset = irtp + 1;
    rtp -= rtp_offset;
    p -=  ip1 * (ip2 + 1) + 1;

    // Function Body 
// ===== 0. INITAIAL SET. 
    vlmin = bzmzro_();
    for (j = 0; j < ncd; ++j)
		tangn[j] = 0.f;
    if (n <= 1)
		return;

    tpc = (set[0] + set[1]) * .5f;
    tlenh = (set[1] - set[0]) * .5f;
    if (tptau <= tpc) {
		i12 = 1;
		r = (tpc - tptau) / tlenh;
    } else {
		i12 = 2;
		r = (tptau - tpc) / tlenh;
    }
    if (itanse[i12] == 1)
		absse = bvabs_(ncd, 1, &tanse[i12*3+1]);
// ====== 1. get approximate tangent from input point seq. 
    np2 = ip;
    if (np2 > n)
		np2 = n;
    kp = k;
    if (kp > np2)
		kp = np2;
    if (ipse == 1) {
		istrt = 1;
		s = tau[1];
    } else {
		istrt = n - np2 + 1;
		s = tau[n];
    }
    for (i = 1; i <= np2; ++i) {
		for (j = 1; j <= ncd; ++j) {
			point[i + j * 6 - 7] = p[(istrt + i - 1 + j * ip2) *  ip1 + 1];
		}
    }
// NP2 IS THE NUMBER OF THE UTILIZED POINT OF P. 
    bkdtkt_(&tau[istrt], np2, kp, work+18);
    iflag=blgint_(&tau[istrt], point, work+18, kp, np2, ncd, ip, ip,work+28, work);
    if(iflag != 1)
		return;

    bleval_(kp, np2, work+18, work, ip, ncd, 1, s, 1, 1, 1, dva);
    rabs = bvabs_(ncd, 1, dva);
	if(rabs > vlmin){
	    for (j = 1; j <= ncd; ++j)
			dva[j - 1] /= rabs;
// NOW APPROXIMATE UNIT TANGENT IN DVA(.) 
// === 2. WHEN TANSE(.,I) GIVEN, DVA AND TANSE ARE MIXED LINEARLY. 
	    if (itanse[i12] == 1) {
			r = r * r * r;
			if (absse >= vlmin) {
				for (j = 1; j <= ncd; ++j) 
					dvb[j - 1] = tanse[j + i12 * 3] / absse;
			    bvivec_(ncd, dva, dvb, r, dva, &t1, &t2);
			}
		}
// ===== 3. IF TAN PLANE GIVEN, OSCULATE TO TAN PLANE. 
	    if (ntp > 0 && (tptau >= ttp[ktp] && tptau <= ttp[ntp + 1])) {
			ble_(ktp, ntp, &ttp[1], &rtp[rtp_offset], irtp, ncd, tptau, 0, dve);
			bvvrt_(dve, dva, dva);
	    }
// NOW 1ST DERIV(UNIT) IN DVA(.) 
// ===== 4. REGULATION OF VECTOR LENGTH 
		if (itanse[i12] == 1)
			rabs += (absse - rabs) * r;
	    for (j = 0; j < ncd; ++j)
			tangn[j] = dva[j] * rabs;
	}else{
	    for (j = 0; j < ncd; ++j)
			tangn[j] = dva[j];
	}
}
