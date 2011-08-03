/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDefint_HH_
#define _MGDefint_HH_

/** @defgroup ALGORITHM (Template) Functions or classes.
 *  @{
 */

#include <math.h>
#include "mg/defintArea.h"

/** 
 *	@brief Integrate f(x)(l=0) or F(x)(l=1, see below) over finite interval (a,b).
 * The DE formula (double exponential formula) is applied.
 * *Note* Original Fortran program codes are from DEFINT of the book 
 * "Fortran77による数値計算ソフトウェア" by M.Sugihara, M.Mori, published by Maruzen K.K.
 *
 *	@retval integral of f(x)(when l=0) or F(x)(when l=1).
 */
template<class func>
double mgDefint(
	func& f,///<Function object for integrand. f(x) returns
		///<the target integration function's value.
	double a,	///<lower bound of integration 
	double b,	///<upper bound of integration
	double eps,	///<absolute error tolerance
	int l=0	///<0 or 1
		///<If l=0 then integrate f(x) over (a,b)
		///<If l=1 then integrate F(x) over (a,b)
		///<In case of l=1 f must be defined as follows:
		///<	IF (-c .LT. y .LT. 0) THEN f(y) = F(b - y) 
		///<	IF ( 0 .LT. y .LE. c) THEN f(y) = F(a - y)
		///<  where c = (b - a) / 2.
){
    if(l!=0 && l!=1)
		return 0.;

	double d;
	int i;
    double wm, wp, fac;

	MGDefintArea::instance().init();//Check if defintArea is initialized.
	int npow=MGDefintArea::instance().m_npow, nend=MGDefintArea::instance().m_nend;
	double *a0=MGDefintArea::instance().m_a0,
		*am=MGDefintArea::instance().m_am,
		*ap=MGDefintArea::instance().m_ap;
	double b0=MGDefintArea::instance().m_b0, *bb=MGDefintArea::instance().m_bb;
	double eps0=MGDefintArea::instance().m_eps0;

	double shfm, shfp;
    fac = (b-a)*.5;
    if(l == 0){
		shfm = (b+a)*.5;
		shfp = shfm;
    }else{
		shfm = 0.;
		shfp = 0.;
    }

    const int lp1 = l + 1;
	const int id=lp1*mgdefintlen-(mgdefintlen+1);
	double epsv;
	eps=fabs(eps);
    if(eps>=eps0)
		epsv = eps;
	else epsv = eps0;
    double epsq = sqrt(epsv) * .2;

    d = a0[l] * fac + shfp;
    double vnew = f(d) * b0;

/*     ---- initial step ---- */
/*          integrate with mesh size = 0.5 and check decay of integrand */
	double h = .5;
	int is=1; is=is<<npow;
    int im = is;
    int km=0, kp=0, nm=0, np=0;
   for(i=is; im<0 ? i>=nend : i<=nend; i+=im){

	if(km<=1){
	    d = am[i+id]*fac + shfm;
	    wm = f(d) * bb[i-1];
	    vnew += wm;
	    if(fabs(wm)<=epsv){
			++km;
			if (km >= 2) nm = i - im;
	    }else
			km = 0;
	}
	if (kp<=1) {
	    d = ap[i+id]*fac + shfp;
	    wp = f(d) * bb[i-1];
	    vnew += wp;
	    if(fabs(wp)<=epsv){
			++kp;
			if (kp>=2) np=i-im;
	    }else
			kp = 0;
	}
	if(km==2 && kp==2)
		break;

    }

    double epsm = 0., epsp = 0.;
    if (nm == 0) {
		nm = nend;
		epsm = sqrt(fabs(wm));
    }

    if(np == 0){
		np = nend;
		epsp = sqrt(fabs(wp));
    }

    if (epsq<epsm)
		epsq = epsm;
    if (epsq<epsp)
		epsq = epsp;

/*     ---- general step ---- */
    double vold = h * fac * vnew;
    for(int mstep = 1; mstep<=npow; ++mstep){
		vnew = 0.;
		int ih = is;
		is /= 2;

		for(i=is; ih<0 ? i>=nm : i<=nm; i+=ih) {
		    d = am[i+id] * fac + shfm;
			vnew += f(d) * bb[i - 1];
		}

		for(i=is; ih<0 ? i>=np : i<=np; i+=ih) {
		    d = ap[i+id] * fac + shfp;
		    vnew += f(d) * bb[i-1];
		}

		vnew = (vold + h*fac*vnew) * .5;
		if((d=vnew-vold,fabs(d))<epsq)
			return vnew;//converged and return

		h *= .5;
		vold = vnew;
    }
    return vnew;
}

/** @} */ // end of ALGORITHM group

#endif
