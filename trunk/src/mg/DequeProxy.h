/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined(__MGDEQUEPROXY_H__)
#define __MGDEQUEPROXY_H__
/** @addtogroup BASE
 *  @{
 */
#include <deque>

template <class T>
class MGDequeProxy{
public:
	typedef typename std::deque<T>::value_type             value_type;
	typedef typename std::deque<T>::reference              reference;
	typedef typename std::deque<T>::const_reference        const_reference;
	typedef typename std::deque<T>::size_type              size_type;
	typedef typename std::deque<T>::difference_type        difference_type;
	typedef typename std::deque<T>::iterator               iterator;
	typedef typename std::deque<T>::const_iterator         const_iterator;
	typedef typename std::deque<T>::reverse_iterator       reverse_iterator;
	typedef typename std::deque<T>::const_reverse_iterator const_reverse_iterator;

	/// Tests two deques for equality.
	inline friend bool operator ==(const MGDequeProxy& deq1, const MGDequeProxy& deq2){
		return *deq1.m_deq == *deq2.m_deq;
	}
	/// Lexicographical comparison.
	inline friend bool operator <(const MGDequeProxy& deq1, const MGDequeProxy& deq2){
		return *deq1.m_deq < *deq2.m_deq;
	}

/// Constructors.
public:
	/// Creates an empty deque.
	MGDequeProxy()
		: m_deq(new std::deque<T>){};

	/// Creates a deque with n elements.
	MGDequeProxy(size_type n)
		: m_deq(new std::deque<T>(n)){};

	/// Creates a deque with n copies of val.
	MGDequeProxy(size_type n, const T& val)
		: m_deq(new std::deque<T>(n, val)){};

	/// The copy constructor.
	MGDequeProxy(const MGDequeProxy& deq)
		: m_deq(new std::deque<T>(*deq.m_deq)){};

	/// Creates a deque with std::deque.
	MGDequeProxy(const std::deque<T>& deq)
		: m_deq(new std::deque<T>(deq)){};

	/// Creates a deque with a copy of a range.
	template <class InputIterator>
	MGDequeProxy(InputIterator first, InputIterator last)
		: m_deq(new std::deque<T>(first, last)){}

/// The destructor.
public:
	~MGDequeProxy(){ delete m_deq;}

/// Operators overload.
public:

	/// The assignment operator.
	MGDequeProxy& operator =(const MGDequeProxy& deq){
		if(this == &deq)
			return *this;
		delete m_deq;
		m_deq = new std::deque<T>(*(deq.m_deq));
		return *this;
	}

	bool operator !=(const MGDequeProxy& deq) const{ return !(*this == deq);}

	/// Returns the n'th element.
	reference operator [](size_type n){ return m_deq->operator[](n);}
	const_reference operator [](size_type n) const{ return m_deq->operator[](n);}

/// Member functions.
public:
	/// Returns an iterator / const_iterator 
	/// pointing to the beginning of the deque.
	iterator begin(){ return m_deq->begin();}
	const_iterator begin() const{ return m_deq->begin();}

	/// Returns an iterator / const_iterator
	/// pointing to the end of the deque.
	iterator end(){ return m_deq->end();}
	const_iterator end() const{ return m_deq->end();}

	/// Returns a reverse_iterator / const_reverse_iterator
	/// pointing to the beginning of the reversed deque.
	reverse_iterator rbegin(){ return m_deq->rbegin();}
	const_reverse_iterator rbegin() const{ return static_cast<const_reverse_iterator>(end());}

	/// Returns a reverse_iterator / const_reverse_iterator
	/// pointing to the end of the reversed deque.
	reverse_iterator rend(){ return m_deq->rend();}
	const_reverse_iterator rend() const{ return static_cast<const_reverse_iterator>(begin());}

	/// Returns the size of the deque.
	size_type size() const{ return m_deq->size();}

	/// Returns the largest possible size of the deque.
	size_type max_size() const{ return m_deq->max_size();}

	/// Returns true if and only if the deque's size is 0.
	bool empty() const{ return m_deq->empty();}

	/// Returns the first element.
	reference front(){ return m_deq->front();}
	const_reference front() const{ return m_deq->front();}

	/// Returns the last element.
	reference back(){ return m_deq->back();}
	const_reference back() const{ return m_deq->back();}

	/// Inserts a new element at the begininng.
	void push_front(const T& val){ m_deq->push_front(val);}
	/// Inserts a new element at the end.
	void push_back(const T& val){ m_deq->push_back(val);}

	/// Removes the first element.
	void pop_front(){ m_deq->pop_front();}
	/// Removes the last element.
	void pop_back(){ m_deq->pop_back();}

	/// Swaps the contents of two deques.
	void swap(MGDequeProxy& deq){ m_deq->swap(*deq.m_deq);}

	/// Inserts val before pos.
	iterator insert(iterator pos, const T& val){
		return m_deq->insert(pos, val);
	}
	/// Inserts the range [first, last) before pos.
	template <class InputIterator>
	void insert(iterator pos, InputIterator first, InputIterator last){
		m_deq->insert(pos, first, last);
	}
	/// Inserts n copies of val before pos.
	void insert(iterator pos, size_type n, const T& val){
		m_deq->insert(pos, n, val);
	}

	/// Erases the element at position pos.
	iterator erase(iterator pos){ return m_deq->erase(pos);}
	/// Erases the range [first, last).
	iterator erase(iterator first, iterator last){ return m_deq->erase(first, last);}

	/// Erases all of the elements.
	void clear(){ m_deq->clear();}

	/// Inserts or erases elements at the end such that the size becomes n.
	void resize(size_type n, const T& val = T()){ m_deq->resize(n, val);}

	/// Returns std::deque object.
	std::deque<T>& std_container(){ return *m_deq;}
	const std::deque<T>& std_container() const{ return *m_deq;}

/// Member data.
private:
	std::deque<T>* m_deq;
};

/** @} */ // end of BASE group
#endif	// __MGDEQUEPROXY_H__