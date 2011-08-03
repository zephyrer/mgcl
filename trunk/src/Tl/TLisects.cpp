/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "Tl/TLisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//mgTLisects is a proprietry class for Face tessellation.
//mgTLisects holds all the intesections of loop(s) and u=const(or v=const)
//line in the parameter space of a face.
//When this holds an intersection with u=const, m_t is u value
//(and t() of a member of mgTLisects is v value), or vise versa.

//Find intersections that are located between (t0,t1) in the m_isects vector.
//Does not include equal.
//When i0>=m_isects.end(), no data between t0 and t1.
//On return, i0=lower_bound() of t0, and i1=upper_bound() of t1.
void mgTLisects::locate(double t0, double t1, iterator& i0, iterator& i1){
	i0=m_isects.begin(); i1=m_isects.end();
	if(i1==i0) return;
	if(m_isects.back()<t0){ i0=i1; return;}
	if(m_isects.front()>t1){ i1=i0; return;}

	//Search t0 and t1's place.
	iterator is=i0, ie=i1;
	i0=std::lower_bound(is,ie,mgTLisect(t0,MGLEPoint(),false));
	i1=std::upper_bound(is,ie,mgTLisect(t1,MGLEPoint(),false));
}

ostream & operator<< (ostream& out, const mgTLisect& is){
	out<<"TLisect::m_t="<<is.m_t<<",m_increace="<<is.m_increase
		<<",m_is=("<<is.m_is<<")";
	return out;
}

ostream & operator<< (ostream& out, const mgTLisects& iss){
	out<<"Tlisects::m_t="<<iss.m_t<<",isect num="<<iss.size()<<endl;
	std::vector<mgTLisect>::const_iterator
		i=iss.m_isects.begin(), ie=iss.m_isects.end();
	size_t j=0;
	for(; i!=ie; i++) out<<j++<<"::"<<(*i)<<endl;
	return out;
}
