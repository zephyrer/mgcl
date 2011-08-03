/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLisect_HH_
#define _mgTLisect_HH_

#include "topo/LEPoint.h"

///mgTLisect is a proprietry class for Face tessellation.
///mgTLisect represents one intesection of Loop and u=const(or v=const)
///line in the parameter space of a face.
///When intersection with u=const, m_t is v value, or vise versa.
class MGCLASS mgTLisect{

public:

MGDECL friend std::ostream & operator<< (std::ostream&, const mgTLisect&);

mgTLisect():m_t(0.),m_increase(0){;};
mgTLisect(double t, const MGLEPoint& is, int increase)
:m_t(t), m_is(is), m_increase(increase){;};

/////////Operator oveload/////////

///Comparison with mgTLisect.
bool operator< (const mgTLisect& is)const{return (m_t<is.m_t);};
bool operator> (const mgTLisect& is)const{return is<*this;};
bool operator<= (const mgTLisect& is)const{return !(*this>is);};
bool operator>= (const mgTLisect& is)const{return !(*this<is);};
bool operator== (const mgTLisect& is)const{return (m_t==is.m_t);};
bool operator!= (const mgTLisect& is)const{return (m_t!=is.m_t);};

///Comparison with double.
bool operator< (double t)const{return (m_t<t);};
bool operator> (double t)const{return m_t>t;};
bool operator<= (double t)const{return m_t<=t;};
bool operator>= (double t)const{return m_t>=t;};
bool operator== (double t)const{return m_t==t;};
bool operator!= (double t)const{return m_t!=t;};

int increase()const{return m_increase;};
const MGLEPoint& isect() const{return  m_is;};
void set_increase(int increase){m_increase=increase;};
double t()const{return m_t;};

private:
	double m_t;	///<u or v value of the rectangle according to the perimeter number.
			///<When perimeter num=0, 2: u value and =1,3:v value.
	MGLEPoint m_is;
	int m_increase;///<Indicates if the loop is increasing or not
			///<at the intersection. =1 if incresing, =-1:decreasing, =0:unknown
};

#endif
