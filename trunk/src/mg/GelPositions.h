/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGelPositions_HH_
#define _MGGelPositions_HH_

class MGGroup;
#include <algorithm>
#include <vector>
#include "mg/GelPosition.h"

class MGPickObjects;

/** @addtogroup GelRelated
 *  @{
 */

///MGGelPosition Container Class.
///MGGelPositions is a class which constains MGGelPosition elements as a vector.
class MGCLASS MGGelPositions{
public:

///	別名定義
	typedef std::vector<MGGelPosition> container_type;
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

	///vector of MGGelPosition.
	container_type m_gelps;

///String output function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGGelPositions&);

////////////Constructor////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGGelPositions(){;}

///Construct MGGelPositions of a MGGelPosition.
MGGelPositions(const MGGelPosition& gelp):m_gelps(1,gelp){;};

///Conversion constructor of MGPickObjects.
MGGelPositions(const MGPickObjects& gelp);

///Copy constructor.
MGGelPositions(const MGGelPositions& obj2);

///Destructor
///~MGGelPositions();

/////////////////operator overloaded////////////////

const_reference operator[](size_type i)const{return m_gelps[i];};
reference operator[](size_type i){return m_gelps[i];};

///Set operation.
MGGelPositions& operator+=(const MGGelPositions& gelps){push_back(gelps);return *this;};
MGGelPositions& operator+=(const MGGelPosition& gelp){push_back(gelp);return *this;};
MGGelPositions& operator-=(const MGGelPositions& gelps){remove(gelps);return *this;};
MGGelPositions& operator-=(const MGGelPosition& gelp){remove(gelp);return *this;};
MGGelPositions& operator-=(const MGAbstractGels& types){remove(types);return *this;};
MGGelPositions& operator&=(const MGGelPositions& gelps){reset_with_common(gelps);return *this;};

////////////Member Function////////////

///Replace this sequence with [first,last).
void assign(const_iterator first, const_iterator last);

/// Return(but does not remove) last element in the group.
/// If list is empty, behavior is undefined.
virtual const_reference back() const{return m_gelps.back();};
virtual reference back(){return m_gelps.back();};

/// Return iterator at the beginning of list.
const_iterator begin() const{return m_gelps.begin();}
iterator begin(){return m_gelps.begin();}

/// clear list, that is, erase all the elements in the MGGelPositions.
void clear(){m_gelps.clear();};

///Return true (1) if there are no items in the MGGelPositions,
/// false(0) otherwise.
bool empty() const{return m_gelps.empty();};

/// Return const_iterator at the end of MGGelPositions.
const_iterator end() const{return m_gelps.end();};
iterator end(){return m_gelps.end();};

/// erase element MGGelPosition gelp. Function's return value is the following iterator
/// of the erased element x.
iterator erase(const MGGelPosition& gelp);

/// erase element x. Function's return value is the following iterator
/// of the erased element x.
iterator erase(iterator x){return m_gelps.erase(x);};

/// erase sequence [first, last). Function's return value is the following iterator
/// of the erased elements.
iterator erase(iterator first, iterator last){return m_gelps.erase(first,last);};

/// erase i-th element.
void erase(size_t i){erase(begin()+i);};

///Find the input MGGelposition.
iterator find(const MGGelPosition& gelp){return std::find(begin(),end(),gelp);};
const_iterator find(const MGGelPosition& gelp)const{return std::find(begin(),end(),gelp);};

///Test if there is a MGGelPosition whose gel is the input gelin
///in this GelPositions' member. If the gel of MGGelPosition is MGShell and gelin
///is MGFace, test is performed if the shell includes the face gelin
///as the shell constituent.
///Returns true if gelin is included in this MGGelPositions.
bool includes(const MGGel* gelin)const;

/// Return(but does not remove) first element in the MGGelPositions.
/// If this vector is empty, behavior is undefined.
virtual const_reference front() const{return m_gelps.front();};
virtual reference front(){return m_gelps.front();};

///insert an element x before the position it.
///Function's return value is the iterator of x after inserted.
iterator insert(iterator it, const MGGelPosition& x){return m_gelps.insert(it,x);};

container_type& object_vector(){return m_gelps;};
const container_type& object_vector()const{return m_gelps;};

/// pop last element.
void pop_back(){m_gelps.pop_back();};

///push elements in gelps at the end. All of the gel pointers are
///transfered to this. On return, gelps will have no gel pointer in it.
void push_back(const MGGelPositions& gelps){
	m_gelps.insert(end(), gelps.begin(), gelps.end());
}
void push_back(const MGGelPosition& gelp){
	m_gelps.push_back(gelp);
}

/// Return reverse_iterator at the beginning of list.
const_reverse_iterator rbegin() const{return m_gelps.rbegin();};
reverse_iterator rbegin(){return m_gelps.rbegin();};

/// Return const_reverse_iterator at the end of list.
const_reverse_iterator rend() const{return m_gelps.rend();};
reverse_iterator rend(){return m_gelps.rend();};

///Resize the gelp.
void resize(size_type n){m_gelps.resize(n);};

///reserve the size n, which are all null.
void reserve(size_t n);

///Remove gelp if found in this.
void remove(const MGGelPosition& gelp){erase(gelp);};

///Remove objects of type from this pickobjects.
void remove(const MGAbstractGels& types);

///Remove gelps from this pickobjects.
void remove(const MGGelPositions& gelps);

///Select objects of specified type from this and reset with them.
void reset_objects(const MGAbstractGels& types);

///replace this with the common objects of this and pobjs2.
void reset_with_common(const MGGelPositions& pobjs2);

///replace this with symmetric_differecne of this and pobj, that is;
///(1) remove the same MGPickObject from this and pobjss.
///(2) append the result pobjs2 to this.
///On return, pobjs2 will have null sequence.
void reset_with_symmetric_difference(MGGelPositions& pobjs2);

///replace the i-th elemnet to pobj.
void reset(size_t i, const MGGelPosition& pobj);

///Select objects of input type from this.
///Function's return value is pickobjects selected.
///This will be unchanged.
MGGelPositions select(const MGAbstractGels& types)const;

///Select the 1st MGCurve from this.
///Function's return value is MGPickObject of MGCurve 1st encountered in this
///MGPickObject sequence. If this did not includes any MGCurve,
///null MGPickOjbect will be returned.
///This will be unchanged.
MGGelPosition select_1st_curve()const;

///Select all the MGCurve from this.
///MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
///in curves.
///This will be unchanged.
void select_curves(MGGelPositions& curves)const;

///Select the 1st MGFSurface from this.
///Function's return value is MGPickObject of MGFSurface 1st encountered in this
///MGPickObject sequence. If this did not includes any MGFSurface,
///null MGPickObject will be returned.
///This will be unchanged.
MGGelPosition select_1st_fsurface()const;

///Select all the MGFSurface from this.
///MGGelPositions of MGFSurface encountered in this MGPickObject sequence will be appended
///in surfaces.
///This will be unchanged.
void select_fsurfaces(MGGelPositions& surfaces)const;

/// Return the number of items that are in the list.
size_type size() const{return m_gelps.size();};

///Test if this is symmetric to gels2.
///Symmetric means:
///(1) number of gels included is the same.
///(2) all of the gels are MGObject and they have the same manifold dimension.
bool symmetric(const MGGelPositions& gels2)const;

};

/** @} */ // end of GelRelated group
#endif
