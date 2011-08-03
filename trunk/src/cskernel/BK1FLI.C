/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//bk1fli_ finds where X is located in xt[.].
//
// *** INPUT * 
//     N,XT(N)....NON-DECREASING REAL NUMBER SEQUENCE OF LENGTH N. 
//     X        A POINT TO BE LOCATED IN XT 
// *** OUTPUT * 
//     Function's return value=THE INDEX LEFT OF XT SUCH THAT ; 
//          XT(LEFT) <= X < XT(LEFT+1)     WHEN 1 <=LEFT < N. 
//          X < XT(1)                      WHEN LEFT=0 
//          XT(N) <= X                      WHEN LEFT=N 
// ******  M E T H O D  ****** 
//  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON 
//  SITUATION THAT IT IS CALLED REPEATEDLY, WITH  X TAKEN FROM 
//  AN INCREASING OR DECREASING SEQUENCE. THIS WILL HAPPEN, E.G., 
//  WHEN A PP FUNCTION IS TO BE GRAPHED. THE FIRST GUESS FOR  LEFT 
//  IS THEREFORE TAKEN TO BE THE VALUE RETURNED AT THE PREVIOUS 
//  CALL AND STORED IN THE  LOCAL VARIABLE  ILO . A FIRST CHECK 
//  ASCERTAINS THAT  ILO .LT. N (THIS IS NECESSARY SINCE THE PRESENT 
//  CALL MAY HAVE NOTHING TO DO WITH THE PREVIOUS CALL). 
//  THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT=ILO  AND 
//  ARE DONE AFTER JUST THREE COMPARISONS. 
//     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE 
//     ISTEP = IHI - ILO 
//  WHILE ALSO MOVING  ILO  AND  IHI  IN THE DIRECTION OF  X , UNTIL 
//                      XT(ILO) .LE. X .LT. XT(IHI) , 
//  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI . 
//  LEFT = ILO  IS THEN RETURNED. 
int bk1fli_(int n, const double *xt, double x){
    // Initialized data 
    static int ilo = 1;
    int istep, middle, ihi;

// Parameter adjustments 
    --xt;

	if(x>=xt[n])
		return n;
	if(x<xt[1])
		return 0;

    ihi = ilo + 1;
	if(ihi>n){
	    ilo = n - 1;
		ihi = n;
	}

    istep=1;
    if(x<xt[ihi]){
	    if(xt[ilo]<=x)
			return ilo;
		//**** NOW X < XT(ILO) . DECREASE  ILO  TO CAPTURE  X . 
		do{
		    ihi = ilo;
			ilo = ihi - istep;
		    istep <<= 1;//Multiply by 2.
		}while(ilo>1 && x<xt[ilo]);
	}else{
		//**** NOW XT(IHI)<=X. INCREASE IHI TO CAPTURE  X . 
		do{
		    ilo = ihi;
			ihi = ilo + istep;
		    istep <<= 1;//Multiply by 2.
		}while(ihi<n && xt[ihi]<=x);
	}
	if(ilo<1)
		ilo=1;
	if(ihi>n)
		ihi=n;

// **** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL. 
	while(1){
	    middle=(ilo+ihi)/2;
		if(middle==ilo)//IT IS ASSUMED THAT MIDDLE=ILO IN CASE IHI=ILO+1. 
			break;
		if(x<xt[middle])
			ihi=middle;
		else
			ilo=middle;
	};
	return ilo;
}

