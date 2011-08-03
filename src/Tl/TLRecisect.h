/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLRecisect_HH_
#define _mgTLRecisect_HH_

#include <vector>
#include "mg/Position.h"
#include "Tl/TLisects.h"

class mgTLRect;

///mgTLRecisect is a proprietry class for Face tessellation.
///mgTLRecisect represents one intesection of Loop and a subdivided (u,v)
///rectangle's perimeter.
class MGCLASS mgTLRecisect{

public:

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLRecisect& ris);

mgTLRecisect():m_perim(-9999){;};
mgTLRecisect(int perim, mgTLisects::iterator is)
:m_perim(perim), m_is(is){;};

/////////Operator oveload/////////
bool operator< (const mgTLRecisect& is)const;
bool operator> (const mgTLRecisect& is)const{return is<*this;}
bool operator<= (const mgTLRecisect& is)const{return !(*this>is);}
bool operator>= (const mgTLRecisect& is)const{return !(*this<is);}
bool operator== (const mgTLRecisect& is)const;
bool operator!= (const mgTLRecisect& is)const{return !operator==(is);}

////////////Member Function/////////////

///Test if this point is going out of or into the rectangle.
bool going_out()const;
bool going_in()const;

mgTLisects::iterator is()const{return m_is;}

///Test if this recisect is in lower range or upper range when mgTLRect is
///subdivided at s.
bool greater_than(
	bool is_u,		///<input if s is u value(true), or not.
	double s);		///<parameter value on a perimeter.
bool less_than(
	bool is_u,		///<input if s is u value(true), or not.
	double s)		///<parameter value on a perimeter.
{return !greater_than(is_u, s);};

int perim() const{return  m_perim;};

///Ask if this isect is processed by createTrimPolygon() or not.
bool processed(){ return m_perim<0;};

///Test if ri2's intersection loop and this are the same.
bool same_loop(const mgTLRecisect& ri2) const;

void set_going_out();
void set_going_in();

///Set the flag that this isect is processed by createTrimPolygon().
void set_processed();

///Return the status of this rec isects.
///=1; going in, =2: going out, =0; unknown;
int status()const;

///Get the (u,v) parameter value of this intersection point.
MGPosition uv(const mgTLRect& rect)const;

private:
	int m_perim;///<Rectanble perimeter number where this intersection lies on.
		///<TLRect::createTrimPolygon() will change the number to -(m_perim+1) using
		///<set_processed() after this rec-isect was processed to generate
		///<the polygon for triangulation.
	mgTLisects::iterator m_is;

friend class mgTLRect;

};

#endif
