/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLLisect_vector_HH_
#define _MGLLisect_vector_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LLisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGLLisect>;
#pragma warning( pop )
#endif

//Forward class declaration.
class MGCurve;
class MGLoop;
class MGInterval;

/// MGLLisect_vector defines a vector of MGLLisect.
/// Used to represent Intersection points of two loops.
class MGCLASS MGLLisect_vector{

public:

typedef std::vector<MGLLisect>::iterator LLiterator;
typedef std::vector<MGLLisect>::const_iterator const_LLiterator;

public:
	std::vector<MGLLisect> m_llisects;///Vector of MGLLisect.

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLLisect_vector& );

////////////Constructor/////////////
MGLLisect_vector();
MGLLisect_vector(const MGLoop& loop);		///Loop

/// Member Function.

/// 交点の全てのコンポーネントを指定して，交点リストに追加
///Add one intersection point to the list.
void append(const MGLLisect& lli);
void append(
	const MGPosition& uv,	///<Parameter (u,v) of the parent face.
	const MGLEPoint& lp1,	///<First loop's point data.
	const MGLEPoint& lp2);	///<Second loop's point data.

/// Adds the MGLLisect_vector to the end of the list.
void append(const MGLLisect_vector& list);

/// Return the number of items that are in the list.
size_t entries() const{return m_llisects.size();};

/// Return first element in the list.
/// If list is empty, behavior is undefined.
const MGLLisect& first() const{return m_llisects.front();};

///Insert MGLLisect at the position i.
void insertAt(LLiterator i, const MGLLisect& llisect)
{m_llisects.insert(i, llisect);};

///Return true if there are no items in the list, false otherwise.
bool empty() const{return m_llisects.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLLisect& last() const{return m_llisects.back();};

private:
	double m_error_square;		///<Error allowed to compute isect and 
								///<end point coincidence.

};

/** @} */ // end of IsectContainer group
#endif
