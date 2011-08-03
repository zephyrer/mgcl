/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLPoint_HH_
#define _MGLPoint_HH_

#include <iosfwd>
#include "mg/MGCL.h"
class MGLEPoint;
//Define MGLPoint Class.

/** @addtogroup TOPORelated
 *  @{
 */

///MGLPoint is to represent Loop's point. This is represented as
///(i, t), where i is the pcell id(i.e. edge number) of in the loop, 
///and t is the parameter value of the curve of the pcell(edge).
class MGCLASS MGLPoint{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLPoint& );

/////////Constructor/////////
MGLPoint():m_i(0), m_t(0.){;};

///Construct from all the necessary data.
MGLPoint(
	size_t i,	///<Loop's edge number.
	double t)	///<Parameter values of i-th pcell curve.
	:m_i(i), m_t(t){;};

///Conversion constructor from MGLEPoint.
MGLPoint(const MGLEPoint& le);

/////////Operator oveload/////////

bool operator< (const MGLPoint& lp)const;
bool operator> (const MGLPoint& lp)const;
bool operator<= (const MGLPoint& lp)const;
bool operator>= (const MGLPoint& lp)const;
bool operator== (const MGLPoint& lp)const;
bool operator!= (const MGLPoint& lp)const{return !operator==(lp);};

/////////Member function/////////

///return loop's edge number.
size_t edge_num()const{return m_i;};

///Return isect data.
double param()const{return m_t;};

private:
	size_t m_i;	///<edge number in the loop.
	double m_t;	///<edge's curve parameter value.

};

/** @} */ // end of TOPORelated group
#endif
