/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// PRIVATE SUB OF BLUMIX, ADD THE FOLLOWING DATA IN THE SEQ , 
//    (KSNEW(NADD),VALNEW(3,NADD),TAUNEW(NADD))  AT VAL(IPOS,3) 
//   ****DATA INSERTED AFTER IPOS**** 
// IFLAG......=0:SUCCESSFULL RETURN, =2:AREA OF VAL EXAUSTED. 
void blumi1_(int *n, int *kseq, double *val, 
	int iv, double *tau, int ipos, int nadd,const int *ksnew,
	const double *valnew, const double *taunew, int *iflag)
{
    int iaft, nmov, i, j, i1, i2;
    // Parameter adjustments 
    --kseq;
    val -= iv + 1;
    --tau;
    --taunew;
    valnew -= 4;
    --ksnew;

    // Function Body 
    if(*n+nadd>iv){
		*iflag=2;
		return;
    }else{
		*iflag=0;
    }
    if(ipos<*n){
		nmov = *n-ipos;
		iaft = *n+nadd;
		for(i=1; i<=nmov; ++i){
		    i1 = iaft-i+1;
			i2 = *n-i+1;
		    kseq[i1] = kseq[i2];
		    for(j=1; j<=3; ++j)
				val[i1+j*iv] = val[i2+j*iv];
		    tau[i1] = tau[i2];
		}
    }
	//   INSERT DATA 
    for(i=1; i<=nadd; ++i){
		i1 = ipos + i;
		kseq[i1] = ksnew[i];
		for(j=1; j<=3; ++j)
			val[i1+j*iv] = valnew[j+i*3];
		tau[i1] = taunew[i];
    }
    *n += nadd;
}

// PRIVATE SUB OF BLUMIX, DELETE NDEL DATA FROM THE SEQ 
void blumi3_(int *n, int *kseq, double *val, 
	int iv, double *tau, int ipos, int ndel)
{
    int iaft, nmov, i, j, i1, i2;
    // Parameter adjustments 
    --kseq;
    val -= iv + 1;;
    --tau;

    nmov = *n-ipos+1;
    iaft = ipos+ndel-1;
    for(i=1; i<=nmov; ++i){
		i1 = iaft+i;
		i2 = ipos+i-1;
		kseq[i2] = kseq[i1];
		for(j=1; j<=3; ++j)
			val[i2+j*iv] = val[i1+j*iv];
		tau[i2] = tau[i1];
    }
    *n -= ndel;
}
