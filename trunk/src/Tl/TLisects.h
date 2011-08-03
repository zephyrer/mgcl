/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLisects_HH_
#define _mgTLisects_HH_

#include <vector>
#include "Tl/TLisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<mgTLisect>;
#pragma warning( pop )
#endif

///mgTLisects is a proprietry class for Face tessellation.
///mgTLisects holds all the intesections of loop(s) and u=const(or v=const)
///line in the parameter space of a face.
///When this holds intersections with u=const, m_t is u value
///(and t() of mgTLisect in mgTLisects is v value), or vise versa.
class MGCLASS mgTLisects{

public:
	
	typedef std::vector<mgTLisect>::iterator iterator;
	std::vector<mgTLisect> m_isects;
		///<intersections of u=m_t(or v=m_t) and the non-perimeter loop
		///<parts of the face's loops. When m_t is u, a member's t() is
		///<v value. Thus these u and v constitute intersection's (u, v) of the face
		///<parameter.
		///<In m_isects the members are sorted in the order of the members' t().

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLisects& is);

/////////////// constructor ////////////////
mgTLisects():m_isects(0){;}
mgTLisects(double t):m_t(t),m_isects(0){;}

/////////Operator oveload/////////
bool operator< (const mgTLisects& is)const{return (m_t<is.m_t);}
bool operator> (const mgTLisects& is)const{return is<*this;}
bool operator<= (const mgTLisects& is)const{return !(*this>is);}
bool operator>= (const mgTLisects& is)const{return !(*this<is);}
bool operator== (const mgTLisects& is)const{return (m_t==is.m_t);}
bool operator!= (const mgTLisects& is)const{return (m_t!=is.m_t);}

///////////// member functions /////////////
mgTLisects::iterator begin(){return m_isects.begin();};
mgTLisects::iterator end(){return m_isects.end();};

///Find intersections that are located between (t0,t1) in the m_isects vector.
///Does not include equal.
///When i0=m_isects.end(), no data between t0 and t1.
///On return, t0<=(i0->t()), and (i1->t())<t1.
void locate(double t0, double t1, iterator& i0, iterator& i1);

void push_back(const mgTLisect& is){m_isects.push_back(is);};
size_t size()const {return m_isects.size();}
void set_t(double t){m_t=t;};
double t()const{return m_t;};

private:
	double m_t;	///<u or v value of the rectangle according to the perimeter number.
		///<When perimeter num=0, 2: u value and =1,3:v value.
};

#endif
