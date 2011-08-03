/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLEPoint_HH_
#define _MGLEPoint_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Unit_vector.h"
#include "topo/Complex.h"

#if defined(MGCL_DLL)
#pragma warning(disable : 4251) /// deriving exported class from non-exported
#endif

class MGLCisect;
class MGLPoint;
class MGLoop;
class MGEdge;

//
//Define MGLEPoint Class.

///MGLEPoint is to represent Loop's point. This is represented as
///(l, i, t), where l is loop pointer, i is the edge's iterator in the loop, 
///and t is the parameter value of the curve of the pcell(edge).
class MGCLASS MGLEPoint{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLEPoint& );

/////////Constructor/////////
MGLEPoint(){;};
MGLEPoint(
	MGComplex::const_pcellItr i,	///<Loop's edge iterator.
	double t)						///<Parameter value of i-th pcell curve.
	:m_i(i), m_t(t){;};

///Construct from the loop, edge number, and the edge's parameter value.
MGLEPoint(
	const MGLoop& lp,
	size_t i,			///<Edge number.
	double t);			///<Parameter value of the edge i.

///Conversion constructor from MGLCisect.
MGLEPoint(const MGLCisect& lci);

///convert loop1's MGLEPoint lep to loop2's MGLEPoint. loop2 must be a copy of loop1.
MGLEPoint(
	const MGLoop& loop1,	///<Loop's edge iterator.
	const MGLEPoint& lep,	///<Parameter value of i-th pcell curve.
	const MGLoop& loop2
);

/////////Operator oveload/////////

///Comparison operator.
bool operator< (const MGLEPoint& lp)const;
bool operator> (const MGLEPoint& lp)const;
bool operator<= (const MGLEPoint& lp)const;
bool operator>= (const MGLEPoint& lp)const;
bool operator== (const MGLEPoint& lp)const;
bool operator!= (const MGLEPoint& lp)const{return !operator==(lp);};

/////////Member function/////////

///Obtain the edge pointer.
const MGEdge* edge() const;
MGEdge* edge_to_update() const;

///return loop's edge number.
size_t edge_num()const;

///Evaluation of the loop at the LEPoint.
///When nderi=0, get the positional data(a parameter (u,v) of the surface)
///at the point.
MGVector eval(size_t nderi=0)const;

///Evaluation of the star curves of the loop at the point t.
///When nderi=0, get a position of the surface at the boundary point.
///The star curve is SurfCurve(face's surface, loop's curve).
///(The star curve has the same world coordinate with the binder curve's, but
///their direction may be opposite. The star curve has always the same direction
///as the loop.)
MGVector eval_star(size_t nderi=0)const;

///Test if two LEPoints are of the same edge.
bool equal_edge(const MGLEPoint& le2) const{return m_i==le2.m_i;};

///Test if two LEPoints are of the same position.
bool equal_position(const MGLEPoint& le2) const;

///Test if this is the same position as P.
bool equal_position(const MGPosition& P) const;

///Compute a vector at the MGLEPoint point that goes inside the face
///and perpendicular to the boundary loop.
MGUnit_vector inner_vector()const;

///test if this is the end point of the loop.
bool is_end_point()const;

///test if this is the start point of the loop.
bool is_start_point()const;

///Return the iterator of the edge.
MGComplex::const_pcellItr iterator()const{return m_i;};

///Get loop pointer.
const MGLoop* loop()const;

///Return edge's curve parameter value.
double param()const{return m_t;};

///Set iteraror of the edge and parameter value of the edge.
void set(
	MGComplex::const_pcellItr i,		///<Iterator of the edge.
	double t);							///<parameter of the edge.

private:
	MGComplex::const_pcellItr m_i;	///<edge iterator in the loop.
	double m_t;						///<edge's curve parameter value.
	

};

///Get the clone of ts.loop()(==te.loop()), and trim th loop from ts to te.
///Function's return value is the trimmed loop.
std::auto_ptr<MGLoop> trim_out_subloop(const MGLEPoint& ts, const MGLEPoint& te);

/** @} */ // end of IsectContainer group
#endif
