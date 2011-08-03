/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "Tl/TLTriangle.h"
#include "Tl/TLFan.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
	
//Test if the edge(i,j) is used or not where i is the index of this edge.
bool mgTLFan::edge_is_used(size_t j)const{
	CEUitr l,ls=m_used_edges.begin(), le=m_used_edges.end();
	l=std::find(ls,le, j);
	return l!=le;
}

//頂点周辺の頂点リストからindexを検索する（後ろから検索）
//If found, iterator of the index be returned.
//If not found , end() will be returned.
mgTLFan::IndexItr mgTLFan::find_aft(size_t index){
	size_t n=size();
	for(int j=n-1; j>=0; j--){
		if(m_indices[j]==index) return begin()+j;
	}
	return end();
}

//Print out indices as "|n0,n1,....
void mgTLFan::print_indices(std::ostream& out,const mgTLTriangle& poly)const{
	size_t n=size();
	size_t nm1=n-1;
	for(size_t i=0; i<n; i++){
		size_t id=(*this)[i];
		out<<id<<"("<<poly[id]<<")";
		if(i<nm1) out<<",";
	}
}

//Set the edge(i,j) as used where i is the index of this fan's vertex.
void mgTLFan::set_edge_used(size_t j){
	EUitr l,ls=m_used_edges.begin(), le=m_used_edges.end();
	l=std::find(ls,le, j);
	if(l==le) m_used_edges.push_back(j);
}

std::ostream& operator<< (std::ostream& out, const mgTLFan& fan){
	size_t n=fan.size();
	out<<"Fan::num of indices="<<n<<"::";
	if(n) out<<fan.m_indices[0];
	for(size_t i=1; i<n; i++) out<<","<<fan.m_indices[i];
	return out;
}
