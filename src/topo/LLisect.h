/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLLisect_HH_
#define _MGLLisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Position.h"
#include "topo/LEPoint.h"

//
//Define MGLLisect Class.

///MGLLisect is to represent two loops intersection point of
///a parent face parameter space.
///Holds two MGLEPoint data of intersection points.
class MGCLASS MGLLisect{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLLisect& );

/////////Constructor/////////

MGLLisect();
MGLLisect(
	const MGPosition& uv,	///<Intersection point data.
	const MGLEPoint& lp1,	///<First loop's LPoint data.
	const MGLEPoint& lp2	///<Second loop's LPoint data.
);

/////////Operator oveload/////////

///Comparison operator.
bool operator< (const MGLLisect& li2)const;
bool operator> (const MGLLisect& li2)const{return li2<(*this);};
bool operator<= (const MGLLisect& li2)const{return !(li2<(*this));};
bool operator>= (const MGLLisect& li2)const{return !((*this)<li2);};
bool operator== (const MGLLisect& li2)const;
bool operator!= (const MGLLisect& li2)const{return !operator==(li2);};

/////////Member function/////////

///Compute distance square of two isect.
double distance_square(const MGLLisect& is2) const;

///Return isect data.
const MGPosition& isect_uv()const{return m_uv;}
const MGLEPoint& isect1()const{return m_is1;};
const MGLEPoint& isect2()const{return m_is2;};

private:
	MGPosition m_uv;	///<parameter value of parent face.
	MGLEPoint m_is1;	///<edge number and param in the first loop.
	MGLEPoint m_is2;	///<edge number and param in the second loop.

};

/** @} */ // end of IsectContainer group
#endif
