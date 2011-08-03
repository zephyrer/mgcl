/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLSPoint_vector_HH_
#define _MGLSPoint_vector_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LSPoint.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGLSPoint>;
#pragma warning( pop )
#endif

/// MGLSPoint_vector defines a vector of MGLSPoint.
/// Used to represent Intersection points of a loop and a surface.
class MGCLASS MGLSPoint_vector{

public:

typedef std::vector<MGLSPoint>::iterator LSiterator;
typedef std::vector<MGLSPoint>::const_iterator const_LSiterator;

public:
	std::vector<MGLSPoint> m_lspoints;///<Vector of MGLSPoint.

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLSPoint_vector& );

////////////Constructor/////////////
MGLSPoint_vector(){;};

/// Member Function.

///Add one intersection point to the list.
bool append(const MGLSPoint& lsp);

/// Adds the MGLSPoint_vector to the end of the list.
void append(const MGLSPoint_vector& vec);

/// Return the number of items that are in the list.
size_t entries() const{return m_lspoints.size();};
size_t size() const{return m_lspoints.size();};

/// Return first element in the list.
/// If list is empty, behavior is undefined.
const MGLSPoint& front() const{return m_lspoints.front();};

///Insert MGLSPoint at the position i.
void insertAt(LSiterator i, const MGLSPoint& llisect)
{m_lspoints.insert(i, llisect);};

///Return true if there are no items in the list, false otherwise.
bool empty() const{return m_lspoints.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLSPoint& back() const{return m_lspoints.back();};

};

/** @} */ // end of IsectContainer group
#endif
