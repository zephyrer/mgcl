/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGAbstractGells_HH_
#define _MGAbstractGells_HH_

#include <vector>
#include "mg/MGCL.h"
#include <iosfwd>
#include "mg/types.h"

/** @addtogroup GelRelated
 *  @{
 */

//
//Define MGAbstractGels Class.

///Is a container of MGAbstractGel, to specify what kind of gels are required.
///MGAbstractGels is a class which constains MGAbstractGel elements as a vector,
///provides OR conditions on the specification of gels.
class MGCLASS MGAbstractGels{
	typedef std::vector<MGAbstractGel> container_type;
public:
//	別名定義
	typedef container_type::reference              reference;
	typedef container_type::const_reference        const_reference;
	typedef container_type::iterator               iterator;
	typedef container_type::const_iterator         const_iterator;
	typedef container_type::size_type              size_type;
	typedef container_type::difference_type        difference_type;
	typedef container_type::value_type             value_type;
	typedef container_type::allocator_type         allocator_type;
	typedef allocator_type::pointer                pointer;
	typedef allocator_type::const_pointer          const_pointer;
	typedef container_type::reverse_iterator       reverse_iterator;
	typedef container_type::const_reverse_iterator const_reverse_iterator;

	///vector of MGAbstractGel.
	std::vector<MGAbstractGel> m_agells;

///String output function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGAbstractGels&);

////////Constructor////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGAbstractGels(){;}

///Construct MGAbstractGels of a MGAbstractGel. This is a conversion contructor.
MGAbstractGels(const MGAbstractGel& agell):m_agells(1,agell){;};

//Copy constructor.
//MGAbstractGels(const MGAbstractGels& obj2);

//Destructor
//~MGAbstractGels();

////////////operator overloaded///////////

///Refer to i-th MGAbstractGel.
const_reference operator[](size_type i)const{return m_agells[i];};
reference operator[](size_type i){return m_agells[i];};

// //////Member Function////////

/// Return(but does not remove) last element in the group.
/// If list is empty, behavior is undefined.
const_reference back() const{return m_agells.back();};
reference back(){return m_agells.back();};

/// Return iterator at the beginning of list.
const_iterator begin() const{return m_agells.begin();}
iterator begin(){return m_agells.begin();}

/// clear list, that is, erase all the elements in the MGAbstractGels.
void clear(){m_agells.clear();};

///Return true (1) if there are no items in the MGAbstractGels,
/// false(0) otherwise.
bool empty() const{return m_agells.empty();};

/// Return const_iterator at the end of MGAbstractGels.
const_iterator end() const{return m_agells.end();};
iterator end(){return m_agells.end();};

/// erase element x.
///Function's return value is the following iterator of the erased element x.
iterator erase(iterator x){return m_agells.erase(x);};

/// erase sequence [first, last).
///Function's return value is the following iterator of the erased elements.
iterator erase(iterator first, iterator last){return m_agells.erase(first,last);};

/// Return(but does not remove) first element in the MGAbstractGels.
/// If this vector is empty, behavior is undefined.
const_reference front() const{return m_agells.front();};
reference front(){return m_agells.front();};

///insert an element x before the position it.
///Function's return value is the iterator of x after inserted.
iterator insert(iterator it, const MGAbstractGel& x){return m_agells.insert(it,x);};

///pop last element.
void pop_back(){m_agells.pop_back();};

///push element x at the end.
void push_back(const MGAbstractGel& x){m_agells.push_back(x);};

///push elements in agells at the end of this.
void push_back(const MGAbstractGels& agells){
	m_agells.insert(end(), agells.begin(), agells.end());
}

///Return reverse_iterator at the beginning of list.
const_reverse_iterator rbegin() const{return m_agells.rbegin();};
reverse_iterator rbegin(){return m_agells.rbegin();};

///Return const_reverse_iterator at the end of list.
const_reverse_iterator rend() const{return m_agells.rend();};
reverse_iterator rend(){return m_agells.rend();};

///Resize the agells.
void resize(size_type n){m_agells.resize(n);};

///Return the number of items that are in the vector.
size_type size() const{return m_agells.size();};

};

/** @} */ // end of GelRelated group
#endif
