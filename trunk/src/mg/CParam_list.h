/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCParam_list_HH_
#define _MGCParam_list_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<double>;
#pragma warning( pop )
#else
#include <list>
#endif

#include <algorithm>
#include "mg/MGCL.h"

//Forward class declaration.
class MGCurve;

/// MGParam_Vector provides a list to store parameters of a curve.
class MGCLASS MGCParam_list{

public:
#if defined(MGCL_DLL)
	typedef MGListProxy<double> container_type;
#else
	typedef std::list<double> container_type;
#endif

container_type m_CPlist;

typedef container_type::iterator Citerator;
typedef container_type::const_iterator const_Citerator;
typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream& operator << (std::ostream&, const MGCParam_list& );

////////////Constructor////////////

/// Constructor of length 0.
explicit MGCParam_list(const MGCurve *curve=0);

//Copy Constructor.
//MGCParam_list(const MGCParam_list& list);

//////////// Destructor.////////////
~MGCParam_list(){;};

//////////// Operator overload.////////////

//Assignment.
//MGCParam_list& MGCParam_list::operator= (const MGCParam_list&);

////////////// Member Function.////////////

/// Adds the parameter to the end of the list.
void append(double param);
void push_back(double param){m_CPlist.push_back(param);};

/// Adds the parameter list to the end of the list.
void append(const MGCParam_list& lst);

///Get the pointer of the first element of the m_CPlist.
Citerator begin(){return m_CPlist.begin();}
const_Citerator begin() const{return m_CPlist.begin();}

///Clear all the elements in m_CPlist.
void clear(){m_CPlist.clear();}

///Returns the pointer to the curve.
const MGCurve* curve() const {return m_curve;}

///Get the pointer of the next of the last element of the m_CPlist.
Citerator end(){return m_CPlist.end();}
const_Citerator end() const{return m_CPlist.end();}

/// Returns the number of items that are in the list.
size_t entries() const{	return size();};
size_t size() const{return m_CPlist.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
Citerator erase(Citerator i){return m_CPlist.erase(i);}

/// Returns(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
double front() const{return m_CPlist.front();};
double first() const{return front();};

///Inserts parameter at the index position i.
///This position must be between zero and the number of items in the list,
/// or behavior is undefined.
void insertAt(Citerator i, double param){m_CPlist.insert(i, param);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
int isEmpty() const{return empty();};
int empty() const{return m_CPlist.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
double last() const{return back();};
double back() const{return m_CPlist.back();};

///Get Lower Bound(this must be sorted)
///if no lower bound return false
bool lower_bound(double param, double& lowerBound)const{
	const_Citerator citer = std::lower_bound(m_CPlist.begin(), m_CPlist.end(), param);
	if(citer == m_CPlist.end())return false;
	lowerBound = *citer;
	return true;
}

///Erase the first element of m_CPlist if not null.
void pop_front(){m_CPlist.pop_front();}

///Erase the last element of m_CPlist if not null.
void pop_back(){m_CPlist.pop_back();}

/// Adds the parameter to the beginning of the list.
void prepend(double param){push_front(param);};
void push_front(double param){m_CPlist.push_front(param);};

///Remove the parameter and return the parameter. If i is not valid, 
/// behavior is undefined.
double removeAt(Citerator i);

///Remove the first parameter int the list and return the parameter.
///If i is not valid, behavior is undefined.
double removeFirst();

///Remove the last parameter in the list and return the parameter.
///If i is not valid, behavior is undefined.
double removeLast();

///Sort the elements in m_CPlist;
void sort(){m_CPlist.sort();};

///erase the same elements as previous element in the sequence.
void unique(){m_CPlist.unique();};

private:
	const MGCurve *m_curve;	///< Curve.
	double m_error;			///<Error in parameter space of the curve.

};

/** @} */ // end of IsectContainer group
#endif
