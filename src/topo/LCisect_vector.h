/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLCisect_vector_HH_
#define _MGLCisect_vector_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LCisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGLCisect>;
#pragma warning( pop )
#endif

//Forward class declaration.
class MGCurve;
class MGLoop;
class MGInterval;

/// MGLCisect_vector defines linked list of MGLCisect.
/// Used to represent Intersection points of a loop and a curve.
class MGCLASS MGLCisect_vector{

public:

typedef std::vector<MGLCisect>::iterator LCiterator;
typedef std::vector<MGLCisect>::const_iterator const_LCiterator;

public:
	std::vector<MGLCisect> m_lcisects;///<Vector of MGLCisect.

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLCisect_vector& );

/// Constructor
MGLCisect_vector();
MGLCisect_vector(const MGLoop& loop);		///Loop

// Operator overload.

const MGLCisect& operator[](size_t i)const{return m_lcisects[i];};
MGLCisect& operator[](size_t i){return m_lcisects[i];};

/// Member Function.

/// Adds the MGLCisect to the end of the list.
void append(const MGLCisect& lcisect);

void append(
	const MGLEPoint& lp,		///<loop's parameter with edge id.
	double t,				///<Curve's parameter value.
	const MGPosition& uv	///<Face's parameter value(u,v) data.
);

/// Adds the MGLCisect_vector to the end of the list.
void append(const MGLCisect_vector& list);

/// Return the number of items that are in the list.
size_t entries() const{return m_lcisects.size();};

/// Return  the first element in the list.
/// If list is empty, behavior is undefined.
const MGLCisect& first() const{return m_lcisects.front();};

///Insert MGLCisect at the position i.
void insertAt(LCiterator i, const MGLCisect& lcisect)
{m_lcisects.insert(i, lcisect);};

///Return true if there are no items in the list, false otherwise.
bool empty() const{return m_lcisects.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLCisect& last() const{return m_lcisects.back();};

///Return the pointer to loop.
const MGLoop* loop() const {return m_loop;};

///Update MGLEPoint in this LCisect_vector.
///This is to update MGLEPoints in this obtained before MGLoop::make_vertex
///and the loop is updated by MGLoop::make_vertex.
///Generally speaking, when make_vertex is invoked after MGLCisect_vector is obtaind,
///the MGLEPoints in the MGLCisect_vector do not contain correct values, since
///a new edge is inserted int the MGComplex's cell vector.
void update_lepoint(
	const MGLEPoint& lep
);

private:
	const MGLoop* m_loop;
	double m_error_square;		///<Error square allowed to compute isect and 
								///<end point coincidence.

};

/** @} */ // end of IsectContainer group
#endif
