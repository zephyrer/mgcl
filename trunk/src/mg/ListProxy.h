/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined(__MGLISTPROXY_H__)
#define __MGLISTPROXY_H__
/** @addtogroup BASE
 *  @{
 */
#include <list>

template <class T>
class MGListProxy{
public:
	typedef typename std::list<T>::value_type             value_type;
	typedef typename std::list<T>::reference              reference;
	typedef typename std::list<T>::const_reference        const_reference;
	typedef typename std::list<T>::size_type              size_type;
	typedef typename std::list<T>::difference_type        difference_type;
	typedef typename std::list<T>::iterator               iterator;
	typedef typename std::list<T>::const_iterator         const_iterator;
	typedef typename std::list<T>::reverse_iterator       reverse_iterator;
	typedef typename std::list<T>::const_reverse_iterator const_reverse_iterator;

	/// Tests two lists for equality. 
	friend bool operator ==(const MGListProxy& ls1, const MGListProxy& ls2){
		return *ls1.m_list == *ls2.m_list;
	}

	/// Lexicographical comparison. 
	friend bool operator <(const MGListProxy& ls1, const MGListProxy& ls2){
		return *ls1.m_list < *ls2.m_list;
	}

	friend void swap(MGListProxy& ls1, MGListProxy& ls2){ ls1.swap(ls2);}

// Constructors.
public:
	/// Creates an empty list.
	MGListProxy()
		: m_list(new std::list<T>){}

	/// The copy constructor.
	MGListProxy(const MGListProxy& ls)
		: m_list(new std::list<T>(*(ls.m_list))){}

	/// Creates a list with n copies of val
	MGListProxy(size_type n, const T& val = T())
		: m_list(new std::list<T>(n, val)){}

	/// Creates a list with std::list.
	MGListProxy(const std::list<T>& ls)
		: m_list(new std::list<T>(ls)){}

/// The destructor.
public:
	~MGListProxy(){ delete m_list;}

/// Operators overload.
public:
	/// The assignment operator.
	MGListProxy& operator =(const MGListProxy& ls){
		if(m_list == ls.m_list)
			return *this;

		delete m_list;
		m_list = new std::list<T>(*(ls.m_list));
		return *this;
	}

/// Member functions.
public:
	/// Returns an iterator pointing to the beginning of the list. 
	iterator begin(){ return m_list->begin();}
	/// Returns a const_iterator pointing to the beginning of the list. 
	const_iterator begin() const{ return m_list->begin();}

	/// Returns an iterator pointing to the end of the list.
	iterator end(){ return m_list->end();}
	/// Returns a const_iterator pointing to the end of the list. 
	const_iterator end() const{ return m_list->end();}

	/// Returns a reverse_iterator pointing to the beginning of the reversed list. 
	reverse_iterator rbegin(){ return m_list->rbegin();}
	/// Returns a reverse_const_iterator pointing to the beginning of the reversed list. 
	const_reverse_iterator rbegin() const{ return static_cast<const_reverse_iterator>(end());}

	/// Returns a reverse_iterator pointing to the end of the reversed list.
	reverse_iterator rend(){ return m_list->rend();}
	/// Returns a reverse_const_iterator pointing to the end of the reversed list. 
	const_reverse_iterator rend() const{ return static_cast<const_reverse_iterator>(begin());}

	/// Returns the first element.
	reference front(){ return m_list->front();}
	const_reference front() const{ return m_list->front();}

	/// Returns the last element.
	reference back(){ return m_list->back();}
	const_reference back() const{ return m_list->back();}

	/// Inserts a new element at the beginning. 
	void push_front(const T& val){ m_list->push_front(val);}
	/// Inserts a new element at the end.
	void push_back(const T& val){ m_list->push_back(val);}

	/// Removes the first element.
	void pop_front(){ m_list->pop_front();}
	/// Removes the last element.
	void pop_back(){ m_list->pop_back();}

	/// Returns the size of the list. 
	size_type size() const{ return m_list->size();}

	/// Returns the largest possible size of the list. 
	size_type max_size() const{ return m_list->max_size();}

	/// Returns true if the list's size is 0.
	bool empty() const{ return m_list->empty();}

	/// Swaps the contents of two lists.
	void swap(MGListProxy& ls){ m_list->swap(*ls.m_list);}

	/// Inserts x before pos.
	iterator insert(iterator pos, const T& val){ return m_list->insert(pos, val);}
	/// Inserts n copies of x before pos.
	void insert(iterator pos, size_type n, const T& val){ m_list->insert(pos, n, val);}

	/// Erases the element at position pos.
	iterator erase(iterator pos){ return m_list->erase(pos);}
	/// Erases the range [first, last).
	iterator erase(iterator first, iterator last){ return m_list->erase(first, last);}

	/// Erases all of the elements.
	void clear(){ m_list->clear();}

	/// Inserts or erases elements at the end such that the size becomes n.
	void resize(size_type n, const T& val = T()){ m_list->resize(n, val);}

	/// pos must be a valid iterator in *this, and ls must be a list that is distinct from *this. 
	/// (That is, it is required that &x != this.) 
	/// All of the elements of ls are inserted before position and removed from pos. 
	/// All iterators remain valid, including iterators that point to elements of ls. 
	/// This function is constant time. 
	void splice(iterator pos, const MGListProxy& ls){ m_list->splice(pos, *ls.m_list);}

	/// pos must be a valid iterator in *this, and i must be a dereferenceable iterator in ls. 
	/// Splice moves the element pointed to by i from ls to *this, 
	/// inserting it before position. 
	/// All iterators remain valid, including iterators that point to elements of ls.
	/// If pos == i or pos == ++i, this function is a null operation. 
	/// This function is constant time. 
	void splice(iterator pos, const MGListProxy& ls, iterator i){ m_list->splice(pos, *ls.m_list, i);}

	/// pos must be a valid iterator in *this, and [first, last) must be a valid range in ls. 
	/// pos may not be an iterator in the range [first, last). 
	/// Splice moves the elements in [first, last) from ls to *this, 
	/// inserting them before position. 
	/// All iterators remain valid, including iterators that point to elements of ls. 
	/// This function is constant time. 
	void splice(iterator pos, const MGListProxy& ls, iterator first, iterator last){
		m_list->splice(pos, *ls.m_list, first, last);
	}

	/// Removes all elements that compare equal to val. 
	/// The relative order of elements that are not removed is unchanged, 
	/// and iterators to elements that are not removed remain valid. 
	/// This function is linear time: 
	/// it performs exactly size() comparisons for equality. 
	void remove(const T& val){ m_list->remove(val);}

	/// Removes all elements *i such that p(*i) is true. 
	/// The relative order of elements that are not removed is unchanged, 
	/// and iterators to elements that are not removed remain valid. 
	/// This function is linear time: it performs exactly size() applications of p. 
	template <class Predicate>
	void remove_if(Predicate p){ m_list->remove_if(p<Predicate>);}

	/// Removes all but the first element in every consecutive group of equal elements. 
	/// The relative order of elements that are not removed is unchanged, 
	/// and iterators to elements that are not removed remain valid. 
	/// This function is linear time: 
	/// it performs exactly size() - 1 comparisons for equality. 
	void unique(){ m_list->unique();}
	/// Removes all but the first element in every consecutive group of equivalent elements, 
	/// where two elements *i and *j are considered equivalent if p(*i, *j) is true. 
	/// The relative order of elements that are not removed is unchanged, 
	/// and iterators to elements that are not removed remain valid. 
	/// This function is linear time: it performs exactly size() - 1 comparisons for equality. 
	template <class BinaryPredicate>
	void unique(BinaryPredicate p){ m_list->unique(p<BinaryPredicate>);}

	/// Both *this and ls must be sorted according to operator<,
	/// and they must be distinct. (That is, it is required that &x != this.) 
	/// This function removes all of ls's elements and inserts them in order into *this. 
	/// The merge is stable; that is, if an element from *this is equivalent to one from ls,
	/// then the element from *this will precede the one from ls. 
	/// All iterators to elements in *this and ls remain valid. 
	/// This function is linear time: 
	/// it performs at most size() + x.size() - 1 comparisons. 
	void merge(MGListProxy& ls){ m_list->merge(*ls.m_list);}

	/// comp must be a comparison function that induces a strict weak ordering 
	/// (as defined in the LessThan Comparable requirements) on objects of type T, 
	/// and both *this and ls must be sorted according to that ordering. 
	/// The lists ls and *this must be distinct. (That is, it is required that &x != this.) 
	/// This function removes all of ls's elements and inserts them in order into *this. 
	/// The merge is stable; that is, if an element from *this is equivalent to one from ls, 
	/// then the element from *this will precede the one from ls. 
	/// All iterators to elements in *this and ls remain valid. 
	/// This function is linear time: it performs at most size() + ls.size() - 1 applications of comp. 
	template <class BinaryPredicate>
	void merge(MGListProxy& ls, BinaryPredicate comp){ m_list->merge(*ls.m_list, comp<BinaryPredicate>);}

	/// Sorts *this according to operator<. 
	/// The sort is stable, that is, the relative order of equivalent elements is preserved. 
	/// All iterators remain valid and continue to point to the same elements. 
	/// The number of comparisons is approximately N log N, where N is the list's size. 
	void sort(){ m_list->sort();}

	/// Reverses the order of elements in the list. 
	/// All iterators remain valid and continue to point to the same elements. 
	/// This function is linear time. 
	void reverse(){ m_list->reverse();}

	/// Returns std::deque object.
	std::list<T>& std_container(){ return *m_list;}
	const std::list<T>& std_container() const{ return *m_list;}

/// Member Data.
private:
	std::list<T>* m_list;
};

/** @} */ // end of BASE group
#endif	// __MGLISTPROXY_H__