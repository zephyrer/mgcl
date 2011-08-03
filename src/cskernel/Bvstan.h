/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
	double *work, double *tangn);
