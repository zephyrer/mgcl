/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLCisect_HH_
#define _MGLCisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Position.h"
#include "topo/LEPoint.h"

//
//Define MGLCisect Class.

class MGLEPoint;
class MGPosition;

///MGLCisect is to represent Loop and curve intersection point of
///a parent face parameter space.
///Holds (lp, t, uv), where lp=loop point, t=curve parameter value, and
///uv=face aprameter value.
class MGCLASS MGLCisect{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLCisect& );

/////////Constructor/////////
MGLCisect();

MGLCisect(
	const MGLEPoint& lp,	///<loop's parameter with edge iterator.
	double t,				///<Curve's parameter value.
	const MGPosition&		///<Face's parameter value(u,v) data.
);

/////////Operator oveload/////////

///Comparison operator.
///THe order is defined as the curve's intersection parameter value ascending
///order if two loops are the same. If two loops are different, loop's address order
///is the order of MGLCisect.
bool operator< (const MGLCisect& li2)const;
bool operator> (const MGLCisect& lci2)const{return lci2<(*this);};
bool operator<= (const MGLCisect& lci2)const{return !(lci2<(*this));};
bool operator>= (const MGLCisect& lci2)const{return !((*this)<lci2);};
bool operator== (const MGLCisect& lci2)const;
bool operator!= (const MGLCisect& lci2)const{return !operator==(lci2);};

/////////Member function/////////

///Compute distance square of two isect.
double distance_square(const MGLCisect& is2) const;

///return loop's edge number.
size_t edge_num()const{return m_lp.edge_num();};

///Get loop pointer.
const MGLoop* loop()const{return m_lp.loop();};

///Get MGLEPoint.
const MGLEPoint& lp() const{return m_lp;};

void set_lepoint(const MGLEPoint& lep){m_lp=lep;};

///Return parameter of the curve.
double t() const {return m_t;}

///Return isect data.
const MGPosition& uv()const{return m_uv;};

private:
	MGLEPoint m_lp;		///<loop's point.
	double m_t;			///<curve's parameter value
	MGPosition m_uv;	///<Face's parameter value

};

/** @} */ // end of IsectContainer group
#endif
