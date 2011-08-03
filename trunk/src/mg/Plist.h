/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPlist_HH_
#define _MGPlist_HH_

// Plist.h
//

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <algorithm>
#include <list>
#include "mg/MGCL.h"

/** @addtogroup BASE
 *  @{
 */

///@brief Defines List of newed object pointers.
///MGPlist is a list of std::auto_ptr<class T>. The member pointers of newed objects will
///be destructed when MGPlist object is destructed. That is, the ownerships of the members
//of MGPlist are transfered to MGPlist, and at the copy and assignment of MGPlist,
///all of the ownerships will be transfered to the copied or assigned new MGPlist.
template <class T>
class MGPlist{
public:
/// Alias.
/// 別名定義 
	typedef typename std::list<T*>::iterator               iterator;
	typedef typename std::list<T*>::const_iterator         const_iterator;
	typedef typename std::list<T*>::reverse_iterator       reverse_iterator;
	typedef typename std::list<T*>::const_reverse_iterator const_reverse_iterator;
	typedef typename std::list<T*>::reference              reference;
	typedef typename std::list<T*>::const_reference        const_reference;
	typedef typename std::list<T*>::size_type              size_type;

/// Member Data
	std::list<T*> m_list;

	///String stream function
	MGDECL friend std::ostream& operator<< (std::ostream& out, const MGPlist<T>& lst);

/// Constructors.

	/// The default constructor.
	MGPlist(){;}

	/// The copy constructor.
	/// The argument rhs will be cleared, since the ownership
	/// of all data of rhs should be transfered to this object.
	///所有権はすべて、新しいオブジェクトに移るので
	///rhsのすべての要素はNULLになります。(sizeは 0となります)
	MGPlist(const MGPlist<T>& rhs) : m_list(rhs.m_list){
		MGPlist<T>& rhsp = const_cast<MGPlist<T>&>(rhs);
		rhsp.m_list.clear();
	}

	/// Create a list whose size is n.
	explicit MGPlist(size_type n) : m_list(n, 0){;}

	/// Create a list with a range. The ownership will be transfered.
	template <class InputIter>
	MGPlist(InputIter first, InputIter last){
		insert(begin(), first, last);
	}

///Destructor
	/// Destroy the object and delete all pointers that this holds.
	~MGPlist();

/// operator overload

	/// assignment operator
	/// The argument rhs will be empty after assignment, 
	/// and transfer the ownership of all pointers to this object.
	/// 代入演算子
	///所有権はすべて、新しいオブジェクトに移るので
	///rhsは空になります)
	MGPlist<T>& operator=(const MGPlist<T>& rhs);

/// member functions

	/// Assignment that takes a range.
	/// The intersection of this object and [first, last) must be empty.
	template <class InputIter>
	void assign(InputIter first, InputIter last){
		iterator cur= begin(), finish = end();
		for(; first != last && cur != finish; ++first, ++cur){
			delete *cur; *cur = *first;
		}
		(first == last) ? erase(cur, finish) : insert(finish, first, last);
	}

	/// Return(but does not remove) last element in the list.
	/// If list is empty, behavior is undefined.
	const_reference back() const{
		return m_list.back();
	}
	reference back(){
		return m_list.back();
	}

	/// Return const_iterator at the beginning of list.
	const_iterator begin() const{
		return m_list.begin();
	}
	iterator begin(){
		return m_list.begin();
	}

	/// clear list, that is, erase all the elements in the list.
	void clear(){
		iterator first = begin(), last = end();
		for(; first != last; ++first)
			delete *first;
		m_list.clear();
	}

	///Return true if there are no items in the list,
	/// false otherwise.
	bool empty() const{
		return m_list.empty();
	}

	/// Return const_iterator at the end of list.
	const_iterator end() const{
		return m_list.end();
	}
	iterator end(){
		return m_list.end();
	}

	/// erase element x.
	iterator erase(iterator x){
		delete *x;
		return m_list.erase(x);
	}
	/// erase sequence [first, last).
	iterator erase(iterator first, iterator last){
		for(; first != last; ++first) delete *first;
		return m_list.erase(first, last);
	}

	///find the position of the object T* found first.
	///If not found, end() will be returned.
	const_iterator find(T* obj) const{
		return std::find(begin(),end(),obj);
	}
	iterator find(T* obj){
		return std::find(begin(),end(),obj);
	}

	/// Return(but does not remove) first element in the list.
	/// If list is empty, behavior is undefined.
	const_reference front() const{
		return m_list.front();
	}
	reference front(){
		return m_list.front();
	}

	///insert an element x before the position it.
	///Function's return value is the iterator of x after inserted.
	iterator insert(iterator it, T* x){
		return m_list.insert(it,x);
	}

	/// Insert given range denoted by [first, last) into one before pos.
	template <class InputIter>
	void insert(iterator pos, InputIter first, InputIter last){
		for(; first != last; ++first) insert(pos, *first);
	}

	///Returns the size of maximum size.
	size_type max_size() const{
		return m_list.max_size();
	}

	/// Equivalent to call std::list::merge().
	void merge(MGPlist<T>& rhs){
		m_list.merge(rhs.m_list);
	}

	/// pop last element.
	void pop_back(){
		delete m_list.back();
		m_list.pop_back();
	}

	/// pop first element.
	void pop_front(){
		delete m_list.front();
		m_list.pop_front();
	}

	/// push element x at the end.
	void push_back(T* x){
		m_list.push_back(x);
	}
	void push_back(MGPlist<T>& x){
		insert(end(), x.begin(), x.end());
		x.m_list.clear();
	}

	/// push element x at the first.
	void push_front(T* x){
		m_list.push_front(x);
	}

	/// Return const_reverse_iterator at the beginning of list.
	const_reverse_iterator rbegin() const{
		return m_list.rbegin();
	}
	reverse_iterator rbegin(){
		return m_list.rbegin();
	}

	/// Return const_reverse_iterator at the end of list.
	const_reverse_iterator rend() const{
		return m_list.rend();
	}
	reverse_iterator rend(){
		return m_list.rend();
	}

	///Release the pointer at th position i.
	///Returned will be the position after the relesed gel.
	iterator release(iterator i){
		return m_list.erase(i);
	}

	/// Release the ownership of elements that point the same as x.
	void remove(T* x){
		m_list.remove(x);
	}

	/// Release the ownership of elements that pred(x) is true, 
	/// where x is an element of the list.
	template <class Pred>
	void remove_if(Pred pred){
		iterator first = begin(), last = end();
		while(first != last){
			iterator next = first;
			++next;
			if(pred(*first)) m_list.erase(first);
			first = next;
		}
	}

	///Remove the T* at the iterator x and return the T*. If x is no valid, 
	/// behavior is undefined.
	///現在の所有権を放棄し、ただのポインタを返す。
	T* removeAt(iterator x){
		T* ret = *x;
		m_list.erase(x);
		return ret;
	}
	
	/// reverse sequence.
	void reverse(){	m_list.reverse();}

	/// Return the number of items that are in the list.
	size_type size() const{ return m_list.size();}

	/// Equivalent to call std::list::sort().
	void sort(){ m_list.sort();}

	void splice(iterator pos, MGPlist<T>& ls){
		m_list.splice(pos, ls.m_list);
	}
	void splice(iterator pos, MGPlist<T>& ls, iterator i){
		m_list.splice(pos, ls.m_list, i);
	}
	void splice(iterator pos, MGPlist<T>& ls, iterator first, iterator last){
		m_list.splice(pos, ls.m_list, first, last);
	}

	/// Equivalent to call std::list::swap().
	void swap(MGPlist<T>& x){
		m_list.swap(x.m_list);
	}

	/// Equivalent to call std::list::unique().
	void unique(){ m_list.unique();}
};

///////////////Implementation////////////

template <class T>
MGPlist<T>::~MGPlist(){
	iterator first = begin(), last = end();
	for(; first != last; ++first) delete *first;
}

template <class T>
MGPlist<T>& MGPlist<T>::operator=(const MGPlist<T>& rhs){
	if(this != &rhs){
		clear();///もとのデータをクリア
		m_list = rhs.m_list;///代入する。
		MGPlist<T>& rhsp=const_cast<MGPlist<T>&>(rhs);
		rhsp.m_list.clear();///代入元をクリア
	}
	return *this;
}

template <class T>
std::ostream& operator<< (std::ostream& out, const MGPlist<T>& lst){
	out << "MGPlist<T>::";
	size_t n=lst.size();
	out<<"number of entries="<<n<<std::endl;
	MGPlist<T>::const_iterator itr; size_t i=0;
	for(itr=lst.begin(); itr!=lst.end(); itr++){
		out<<i++<<":"<<(*itr)<<std::endl;
		if(*itr) out << (**itr) << std::endl;
	}
	return out;
}

/** @} */ // end of BASE group
#endif
