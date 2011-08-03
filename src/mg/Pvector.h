/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPvector_HH_
#define _MGPvector_HH_

#if _MSC_VER > 1000
#pragma once
#endif /// _MSC_VER > 1000

#include <algorithm>
#include <vector>
#include "mg/MGCL.h"
/** @addtogroup BASE
 *  @{
 */

///@brief Defines Vector of newed object pointers.
///MGPvector is a list of std::auto_ptr<class T>. The member pointers of newed objects will
///be destructed when MGPList object is destructed. That is, the ownerships of the members
//of MGPvector are transfered to MGPvector, and at the copy and assignment of MGPvector,
///all of the ownerships will be transfered to the copied or assigned new MGPvector.
template <class T>
class MGPvector{
public:
/// Alias.
///	別名定義
	typedef typename std::vector<T*>::iterator               iterator;
	typedef typename std::vector<T*>::const_iterator         const_iterator;
	typedef typename std::vector<T*>::reverse_iterator       reverse_iterator;
	typedef typename std::vector<T*>::const_reverse_iterator const_reverse_iterator;
	typedef typename std::vector<T*>::reference              reference;
	typedef typename std::vector<T*>::const_reference        const_reference;
	typedef typename std::vector<T*>::size_type              size_type;
	
///Member data
	std::vector<T *> m_vector;

	/// Debug out
	template <class T>
	MGDECL friend std::ostream& operator<< (std::ostream& out, const MGPvector<T>& vector);

//////////// Constructor. /////////

	///void constructor
	MGPvector(){;}

	///default constructor
	/// Create a vector whose size is n.
	///長さlenのベクトルを作成します。
	///すべての要素はNULLで初期化されます。
	explicit MGPvector(size_type len) : m_vector(len, 0){;}

	///Construct MGPvector of length 1 that contains ptr
	MGPvector(T* ptr) : m_vector(1, ptr){;}

	///copy constructor
	/// The argument rhs will be cleared, since the ownership
	/// of all data of rhs should be transfered to this object.
	///所有権はすべて、新しいオブジェクトに移るので
	///rhsのすべての要素はNULLになります。(sizeは変わりません)
	MGPvector(const MGPvector<T>& rhs)
		:m_vector(rhs.m_vector){
		MGPvector<T>& rhsp=const_cast<MGPvector<T>&>(rhs);
		rhsp.m_vector.clear();
	}

	///Conversion constructor.
	MGPvector(const std::vector<const T*>& rhs)
	:m_vector(rhs.size()){
		size_t n=rhs.size();
		for(size_t i=0; i<n; i++){
			m_vector[i]=new T(*(rhs[i]));
		}		
	}

	/// Create a vector with given range [first, last).
	template <class InputIter>
	MGPvector(InputIter first, InputIter last){
		insert(end(), first, last);
	}

//////////// destructor ////////////

	/// Destroy the object and delete all pointers that this object holds.
	~MGPvector();

////////////  operator overload ////////////

	/// Subscript access.
	const_reference operator[](size_type pos) const{
///		if(pos>=size()) return 0;
		return m_vector[pos];
	}
	reference operator[](size_type pos){
///		if(pos>=size()) return 0;
		return m_vector[pos];
	}

	/// assignment operator
	/// The argument rhs will be empty after assignment, 
	/// and transfer the ownership of all pointers to this object.
	///代入演算子
	///所有権はすべて、新しいオブジェクトに移るので
	///rhsのすべての要素はNULLになります。(sizeは変わりません)
	MGPvector<T>& operator=(const MGPvector<T>& rhs);

//////////// Member function ////////////

	/// Assignment that takes a range.
	/// The intersection of this object and [first, last) must be empty.
	template <class InputIter>
	void assign(InputIter first, InputIter last){
		iterator cur= begin();
		for(; first != last && cur != end(); ++first, ++cur){
			delete *cur; *cur = *first;
		}
		(first == last) ? erase(cur, end()) : insert(end(), first, last);
	}

	void assign(size_type pos, T* ptr){
		delete m_vector[pos];
		m_vector[pos]=ptr;
	}

	/// Equivalent to call vector::at(n).
	const_reference at(size_type n) const{
		return m_vector.at(n);
	}
	reference at(size_type n){
		return m_vector.at(n);
	}

	/// Return the reference to the last element in the sequence.
	/// If the vector is empty, behavior is undefined.
	const_reference back() const{
		return m_vector.back();
	}
	reference back(){
		return m_vector.back();
	}

	/// Return const_iterator that points to the first element in the sequence.
	const_iterator begin() const{
		return m_vector.begin();
	}
	/// Return iterator that points to the first element in the sequence.
	iterator begin(){
		return m_vector.begin();
	}

	/// Equivalent to call vector::capacity().
	size_type capacity() const{
		return m_vector.capacity();
	}

	/// clear Sequence, that is, erase all the elements in the sequence.
	void clear(){
		for(iterator i = begin(); i != end(); ++i)
			delete *i;
		m_vector.clear();
	}

	///Return true (1) if there are no items in the vector,
	/// false(0) otherwise.
	bool empty() const{
		return m_vector.empty();
	}

	/// Return const_iterator that points one past the last element.
	const_iterator end() const{
		return m_vector.end();
	}
	/// Return iterator that points one past the last element.
	iterator end(){
		return m_vector.end();
	}

	///erase the element at x.
	iterator erase(iterator x){
		delete *x;
		return m_vector.erase(x);
	}

	///erase i-th element x.
	iterator erase(size_type i){
		iterator del = m_vector.begin() + i;
		delete *del;
		return m_vector.erase(del);
	}

	/// erase sequence [first, last).
	iterator erase(iterator first, iterator last){
		iterator cur = first;
		for(; cur != last; ++cur) delete *cur;
		return m_vector.erase(first, last);
	}

	/// Return the reference to first element in the vector.
	/// If the vector is empty, behavior is undefined.
	const_reference front() const{
		return m_vector.front();
	}
	reference front(){
		return m_vector.front();
	}

	///ただのポインタを返す。所有権は放棄しないので
	///このポインタを解放してはならない。
	const_reference get(size_type i) const{
		return this->operator[](i);
	}
	reference get(size_type i){
		return this->operator[](i);
	}

	/// Insert ptr into vector at pos.
	iterator insert(iterator pos, T* ptr){
		return m_vector.insert(pos, ptr);
	}
	/// Insert the range [first, last) into the vector at pos.
	template <class InputIter>
	void insert(iterator pos, InputIter first, InputIter last){
		for(; first != last; ++first){
			pos = insert(pos, *first);
			++pos;
		}
	}

	///Test if i-th element is null.
	bool is_null(size_type i)const{
		return m_vector[i]==0;
	}

	///Returns the size of maximum size.
	size_type max_size() const{
		return m_vector.max_size();
	}

	/// pop last element.
	void pop_back(){
		delete m_vector.back();
		m_vector.pop_back();
	}

	/// push element x at the end.
	void push_back(T* x){
		m_vector.push_back(x);
	}

	/// push element x at the end.
	///All the ownership of the elements in dst are transfered
	///to this Pvector.
	void push_back(MGPvector<T>& dst){
		m_vector.insert(m_vector.end(), dst.m_vector.begin(), dst.m_vector.end());
		dst.m_vector.clear();
	}

	/// Return const_reverse_iterator that points the end of the sequence.
	const_reverse_iterator rbegin() const{
		return m_vector.rbegin();
	}
	/// Return reverse_iterator that points the end of the sequence.
	reverse_iterator rbegin(){
		return m_vector.rbegin();
	}

	/// Return const_reverse_iterator that points to one before the
	/// first element in the vector.
	const_reverse_iterator rend() const{
		return m_vector.rend();
	}
	/// Return reverse_iterator that points to one before the
	/// first element in the vector.
	reverse_iterator rend(){
		return m_vector.rend();
	}

	///現在の所有権を放棄し、ただのポインタを返す。
	///この処理の後はget(i)==nullとなる。
	///ベクトルの要素自体は削除されません。
	T* release(size_type i){
		T* ret = m_vector[i];
		m_vector[i] = 0;
		return ret;
	}
	///現在の所有権をすべて放棄する。
	///この処理の後はget(i)==nullとなる。
	///ベクトルの要素自体は削除されない。
	void release_all(){
		m_vector.assign(size(), static_cast<T*>(0));
	}

	///Remove the T* and return the T*. If i is no valid, 
	/// behavior is undefined.
	T* removeAt(iterator i){
		T* ret = *i;
		m_vector.erase(i);
		return ret;
	}
	///Remove the i-th T* in the vector and return the T*.
	///If i is not valid, the behavior is undefined.
	T* removeAt(size_type i){
		T* ret = m_vector[i];
		iterator del = m_vector.begin() + i;	
		m_vector.erase(del);
		return ret;
	}

	///reverse the sequence of the elements.
	void reverse_sequence(){
		size_type i,n=size();
		size_type nhalf=n/2;
		for(i=0; i<nhalf; i++){
			size_type nmim1=n-i-1;
			T* pointer=m_vector[nmim1];
			m_vector[nmim1]=m_vector[i];
			m_vector[i]=pointer;
		}
	}

	/// Preallocate memory for specified number of elements
	/// if necessary.
	void reserve(size_type n){
		m_vector.reserve(n);
	}

	/// delete the i-th element and replace it with ptr.
	void reset(size_type i, T* ptr = 0){
		delete m_vector[i];
		m_vector[i] = ptr;
	}

	///Resize the vector.
	///When size is enlarged, enlarged part will contain null pointer.
	void resize(size_type n);

	/// Return the number of items that are in the sequence.
	size_type size() const{
		return m_vector.size();
	}

	void swap(MGPvector<T>& x){
		m_vector.swap(x.m_vector);
	}
};

///////////////////////////////////////////////////////////////////////////////
/// Implementation. 
///////////////////////////////////////////////////////////////////////////////

template <class T>
MGPvector<T>::~MGPvector(){
	for(iterator i = begin(); i != end(); ++i) delete *i;
}

template <class T>
MGPvector<T>& MGPvector<T>::operator=(const MGPvector<T>& rhs){
	if(this != &rhs){
		clear(); ///もとのデータをクリア
		m_vector = rhs.m_vector; ///代入する。
		MGPvector<T>& rhsp=const_cast<MGPvector<T>&>(rhs);
		rhsp.m_vector.clear(); ///代入元をクリア
	}
	return *this;
}

///Resize the vector.
///When size is enlarged, enlarged part will contain null pointer.
template <class T>
void MGPvector<T>::resize(size_type n){
	size_t oldn=m_vector.size();
	if(oldn<n) m_vector.resize(n,0);
	else if(oldn>n){
		for(size_t i=n; i<oldn; i++) delete m_vector[i];
		m_vector.resize(n);
	}
}

template <class T>
std::ostream& operator<< (std::ostream& out, const MGPvector<T>& vector){
	out << "MGPvector<T>::";
	size_t n=vector.size();
	out<<"number of entries="<<n<<std::endl;
	MGPvector<T>::const_iterator itr; size_t i=0;
	for(itr=vector.begin(); itr!=vector.end(); itr++){
		out<<i++<<":"<<(*itr)<<":";
		if(*itr) out << (**itr) << std::endl;
	}
	return out;
}

/** @} */ // end of BASE group
#endif
