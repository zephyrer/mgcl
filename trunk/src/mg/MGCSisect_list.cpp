/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/CSisect_list.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCSisect_list defines linked list of MGCSisect.
// Used to represent Intersection points of two curves.

#define ERROR_EXPAND 5.
// Constructor
MGCSisect_list::MGCSisect_list(const MGCurve* crv, const MGFSurface* srf)
:m_curve(crv), m_surface(srf)
{
	m_errort=m_erroru=m_errorv=0.;
	if(crv) m_errort=crv->param_error()*ERROR_EXPAND;
	if(srf){
		m_erroru=srf->param_error_u()*ERROR_EXPAND;
		m_errorv=srf->param_error_u()*ERROR_EXPAND;
	}
}

//Copy Constructor.

// Destructor.

// Operator overload.

//Assignment.

// Member Function.

void MGCSisect_list::append(const MGCSisect& isect){
// Adds the MGCSisect to the end of the list.
//End points will be prefered.
	double t0=m_curve->param_s(), t1=m_curve->param_e(), ter;
	double uer,ver;

	CSiterator itr;
	for(itr=begin(); itr!=end(); itr++){
		isect.distance(*itr,ter,uer,ver);
		if(ter<=m_errort && uer<=m_erroru && ver<=m_errorv){
			double s=(*itr).param_curve(), t=isect.param_curve();
			double slen, tlen;
			if(s<=t){
				slen=s-t0; tlen=t1-t;
			}else{
				slen=t1-s; tlen=t-t0;
			}
			if(slen<=tlen) return;
			else itr=erase(itr);
		}
	}

	push_back(isect);
}

// 全てのコンポーネントを指定して交点を生成する
void MGCSisect_list::append(
		const MGPosition& point,		//intersection point.
		double t,				//Curve's parameter value.
        const MGPosition& uv,	//Surface's parameter values.
		const MGCSRELATION rl	//Curve and Surface relation
		)
{
	append(MGCSisect(point,t,uv,rl));
}

void MGCSisect_list::append(const MGCSisect_list& list){
// Adds the MGCSisect_list to the end of the list.
	const_CSiterator i;
	for(i=list.begin(); i!=list.end(); i++) append(*i);
}

MGCSisect MGCSisect_list::removeAt(CSiterator i){
//Remove the MGCSisect and return the MGCSisect. If i is no valid, 
// behavior is undefined.
	MGCSisect isect=*i;
	erase(i);
	return isect;
}

MGCSisect MGCSisect_list::removeFirst(){
//Remove the first MGCSisect int the list and return the MGCSisect.
//If i is not valid, behavior is undefined.
	MGCSisect isect=front();
	pop_front();
	return isect;
}

MGCSisect MGCSisect_list::removeLast(){
//Remove the first MGCSisect int the list and return the MGCSisect.
//If i is not valid, behavior is undefined.
	MGCSisect isect=back();
	pop_back();
	return isect;
}
