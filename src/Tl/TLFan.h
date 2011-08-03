/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLFan_HH_
#define _mgTLFan_HH_

#include <algorithm>
#include <vector>
#include <iosfwd>

#include "mg/MGCL.h"
#if defined(MGCL_DLL)
#include "mg/DequeProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGDequeProxy<size_t>;
#pragma warning( pop )
#else
#include <deque>
#endif

class mgTLTriangle;

///////////// mgTLFan /////////////

// private class for tessellation.

///mgTLFan is a point list to constitue a fan.
///It is always a member of mgTLFans that is temporary data to generate
///mgTLTriangles, and contains id of mgTLPoints that is a vector of MGPosition uv
///of a surface parameter. Both of mgTLPoints and mgTLTriangles are members of
///mgTLData.
///mgTLFan does not include its start point id since the id in mgTLFans is the 1st id.
///Let mgTLFan& fani=mgTLFans[i], then (i,fani[0], fani[1], ..., fani[n-1]) constitutes
///a fan(all of them are id of mgTLPoints). The actual triangles are:
///(i,fani[0], fani[1]), (i,fani[1], fani[2]), ..., (i,fani[i], fani[i+1]), ...,
///(i, fani[n-2], fnai[n-1]).
///See Les Piegel and Wayne Tiller's paper about this point list.
///"Geometry-based triangulation of trimmed NURBS surfaces", Computer-Aided Desigh,
///Vol.30, No.1, pp.11-16, 1998.
class MGCLASS mgTLFan{
private:
#if defined(MGCL_DLL)
	typedef MGDequeProxy<size_t> mgTLdeqIndex;
#else
	typedef std::deque<size_t> mgTLdeqIndex;
#endif

	bool m_used;	///<flag to indicate if this vertex is used or not.
	mgTLdeqIndex	m_indices;///<includes vertices from the second points to the last.
	std::vector<size_t> m_used_edges;
					///<Include j's array if the edge(i,j) is a used edge in mgTLFans.
					///<j is always greater than i, and the order of j's is undefined.
					///<Here i is this fan's index.				

public:
	typedef mgTLdeqIndex::iterator IndexItr;
	typedef mgTLdeqIndex::const_iterator CIndexItr;
	typedef mgTLdeqIndex::reverse_iterator ritr;
	typedef std::vector<size_t>::iterator EUitr;///<Edge Used iterator.
	typedef std::vector<size_t>::const_iterator CEUitr;///<Edge Used iterator.

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLFan& fan);

///////////////Constructor///////////////

	mgTLFan():m_used(false){;};
	mgTLFan(size_t v1):m_indices(1,v1){;};
	mgTLFan(size_t v1, size_t v2):m_used(false),m_indices(2){m_indices[0]=v1;m_indices[1]=v2;};
	mgTLFan(mgTLdeqIndex& index):m_used(false),m_indices(index){;};

//////////// Operator overload ///////////////
	size_t operator[](size_t i)const{return m_indices[i];};

/////////////member functions////////////

	size_t back()const{return m_indices.back();};
	IndexItr begin(){return m_indices.begin();};
	CIndexItr begin()const{return m_indices.begin();};
	
	///Test if the edge(i,j) is used or not where i is the index of this edge.
	bool edge_is_used(size_t j)const;

	IndexItr end(){return m_indices.end();};
	CIndexItr end()const{return m_indices.end();};
	void erase(IndexItr iter){m_indices.erase(iter);};

	///頂点周辺の頂点リストからindexを検索する（前から検索）
	CIndexItr find(size_t index)const{
		return std::find(m_indices.begin(), m_indices.end(), index);
	}
	///頂点周辺の頂点リストからindexを検索する（前から検索）
	IndexItr find(size_t index){
		return std::find(m_indices.begin(), m_indices.end(), index);
	}

	///頂点周辺の頂点リストからindexを検索する（後ろから検索）
	IndexItr find_aft(size_t index);

	size_t front()const{return m_indices.front();};

	///Insert the index before the position iter.
	IndexItr insert(IndexItr iter, size_t index){
		return m_indices.insert(iter, index);
	};

	const mgTLdeqIndex& indices()const{return m_indices;};
	void push_back(size_t index){m_indices.push_back(index);};
	void push_front(size_t index){m_indices.push_front(index);};
	void pop_back(){m_indices.pop_back();};
	void pop_front(){m_indices.pop_front();};

	///Print out indices as "|n0,n1,....
	void print_indices(std::ostream& out,const mgTLTriangle& poly)const;

	ritr rbegin(){return m_indices.rbegin();};
	ritr rend(){return m_indices.rend();};

	///Set this vertex as used.
	void set_vertex_used(){m_used=true;};

	///Set the edge(i,j) as used where i is the index of this fan's vertex.
	void set_edge_used(size_t j);

	size_t size()const{return m_indices.size();};

	bool vertex_is_used()const{return m_used;};

};

#endif
