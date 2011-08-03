/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
	int *kplane, double *g, double *cntr, double *devmax);
