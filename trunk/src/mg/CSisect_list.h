/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCSisect_list_HH_
#define _MGCSisect_list_HH_

/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/CSisect.h"

#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<MGCSisect>;
#pragma warning( pop )
#else
#include <list>
#endif

//Forward class declaration.
class MGCurve;
class MGFSurface;

/// MGCSisect_list defines linked list of MGCSisect.
/// Used to represent Intersection points of a curve and a surface.
class MGCLASS MGCSisect_list{

public:
#if defined(MGCL_DLL)
	typedef MGListProxy<MGCSisect> container_type;
#else
	typedef std::list<MGCSisect> container_type;
#endif

container_type m_CSilist;

typedef container_type::iterator CSiterator;
typedef container_type::const_iterator const_CSiterator;

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGCSisect_list& );

////////////// Constructor////////////
explicit MGCSisect_list(const MGCurve *crv=NULL, const MGFSurface *srf=NULL);

/// Destructor.
~MGCSisect_list(){;};

////////////// Member Function. ////////////

void append(const MGCSisect& isect);
void push_back(const MGCSisect& isect){m_CSilist.push_back(isect);};

/// 全てのコンポーネントを指定して交点を生成する
void append(
		const MGPosition& point,		///<intersection point.
		double t,				///<Curve's parameter value.
        const MGPosition& uv,	///<Surface's parameter values.
		const MGCSRELATION rl=MGCSREL_UNKNOWN
								///<Curve and Surface relation
	);

void append(const MGCSisect_list& list);

///Get the pointer of the first element of the m_CSilist.
CSiterator begin(){return m_CSilist.begin();}
const_CSiterator begin() const{return m_CSilist.begin();}

///Clear all the elements in m_CSilist.
void clear(){m_CSilist.clear();}

///Return the pointer to the curve.
const MGCurve* curve() const {return m_curve;};

///Get the pointer of the next of the last element of the m_CSilist.
CSiterator end(){return m_CSilist.end();}
const_CSiterator end() const{return m_CSilist.end();}

/// Return the number of items that are in the list.
size_t entries() const{	return size();};
size_t size() const{return m_CSilist.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
CSiterator erase(CSiterator i){return m_CSilist.erase(i);}

/// Return(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
const MGCSisect& first() const{return front();};
const MGCSisect& front() const{return m_CSilist.front();};

///Insert MGCSisect at the iterator i.
void insertAt(CSiterator i, const MGCSisect& isect)
{m_CSilist.insert(i, isect);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
int isEmpty() const{return empty();};
int empty() const{return m_CSilist.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGCSisect& last() const{return m_CSilist.back();};
const MGCSisect& back() const{return m_CSilist.back();};

///Erase the last element of m_CSilist if not null.
void pop_back(){m_CSilist.pop_back();}

///Erase the first element of m_CSilist if not null.
void pop_front(){m_CSilist.pop_front();}

/// Adds the parameter to the beginning of the list.
void prepend(const MGCSisect& param){push_front(param);};
void push_front(const MGCSisect& isect){m_CSilist.push_front(isect);};

///Remove the parameter and return the parameter. If i is not valid, 
/// behavior is undefined.
MGCSisect removeAt(CSiterator i);

///Remove the first MGCSisect int the list and return the MGCSisect.
///If i is not valid, behavior is undefined.
MGCSisect removeFirst();

///Remove the last MGCSisect in the list and return the MGCSisect.
///If i is not valid, behavior is undefined.
MGCSisect removeLast();

///Return the pointer to the surface.
const MGFSurface* surface() const {return m_surface;};

private:
	const MGCurve *m_curve;		///< Curve.
	const MGFSurface *m_surface;	///< Surface.
	double m_errort;///<error to regard same curve point in parameter space.
	double m_erroru;///<error to regard same surface point in u-parameter space.
	double m_errorv;///<error to regard same surface point in v-parameter space.

};

/** @} */ // end of IsectContainer group
#endif
