/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLisetsList_HH_
#define _mgTLisetsList_HH_

#include "Tl/TLisects.h"
#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<mgTLisects>;
#pragma warning( pop )
#else
#include <list>
#endif

class mgTLparameter;

///mgTLisectsList is a proprietry class for Face tessellation.
///mgTLisectsList holds all the necessary mgTLisects along u or v 
///in the parameter space of a face.
///This class is to eliminate the same value's isect_1D computation.
///Once isect_1D is invoked for some u or v value, the intersections are
///hold in this class.
class MGCLASS mgTLisectsList{

public:
#if defined(MGCL_DLL)
	typedef MGListProxy<mgTLisects> container_type;
#else
	typedef std::list<mgTLisects> container_type;
#endif

	typedef container_type::iterator issListItr;
	typedef container_type::const_iterator CissListItr;
	container_type m_isectsList;
	///<intersections of u=m_t(or v=m_t) and the non-perimeter loop
	///<parts of the face's loops.

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLisectsList& is);

///////////// Constructor ////////////
mgTLisectsList():m_isectsList(0){;};

///////////// Member function /////////////

///Search mgTLisects whose t()==tin. If not found compute intersections
///and add one member of mgTLisects whose t() is tin.
///Function's return value is the iterator of the mgTLisects in mgTLisectsList.
issListItr get(
	mgTLparameter& param,	///<tessellation parameter.
	double tin,	///<face's parameter u or v according to kcod.
	size_t kcod	///<Coordinate kind, =0: tin is u value, =1: tin is v value.
);

void clear(){m_isectsList.clear();};
double middle()const{return m_middle;};
double minimum()const{return m_isectsList.begin()->t();};
double maximum()const{return m_isectsList.rbegin()->t();};

void push_back(const mgTLisects& iss){m_isectsList.push_back(iss);};
void push_front(const mgTLisects& iss){m_isectsList.push_front(iss);};
void push_back(double t){m_isectsList.push_back(mgTLisects(t));};
void push_front(double t){m_isectsList.push_front(mgTLisects(t));};
void set_middle(double t){ m_middle=t;};

private:
	double m_middle; ///<(minimum()+maximum())*0.5 will be stored in TLinit.
};

///Compute intersection points of 1D sub curves of the original loop.
///Parameter values of intersection points(MGLEPoints) will be returned.
///This is for tessellation and intersection with perimeter boudary edges will
///be excluded.
MGDECL void mgTLisect1D(
	mgTLparameter& param,	///<tessellation parameter.
	const MGLoop& lp,	///<Target Loop.
	const std::vector<bool>& edgPerimj,///<lp's edgPerim[j].
		///<if lp is an inner boundary, this is not used.
		///<Used only when lp is outer boundary or perimeter boundary.
	double f,		///<Coordinate value
	size_t kcod,	///<Coordinate kind of the data f, u(kcod=0) or v(=1).
	mgTLisects& iss	///<mgTLisect vector will be output.
					///<Must be initialized before use.
);

#endif
