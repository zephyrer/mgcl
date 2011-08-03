/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD504.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/BVertex.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD504.h"
#include "mgiges/Iges504EdgeListMap.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD504 is the class for Iges parameter data type 504(EDGE list).

// Constructors.

//! Constructs an object of class MGIgesPD504.
MGIgesPD504::MGIgesPD504(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(EDGE,DEpointer),m_edges(1,MGIges504Edge()){
}

//Read in parameter data from string stream data.
void MGIgesPD504::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_edges;
	get_integer(pDelimeter,pdstream,num_edges);
	m_edges.resize(num_edges+1);
	for(int i=1; i<=num_edges; i++){
		int curveDE;
		get_DEpointer(pDelimeter,pdstream,curveDE);
		int Svertex_list;//pointer to the DE of the VERTEX List Entry(MGPD502) for the start vertex.
		get_DEpointer(pDelimeter,pdstream,Svertex_list);
		int Svertex;//List index of the start vertex in m_Svertex_list DE.
		get_integer(pDelimeter,pdstream,Svertex);
		int Tvertex_list;//pointer to the DE of the VERTEX List Entry(MGPD502) for the terminate vertex.
		get_DEpointer(pDelimeter,pdstream,Tvertex_list);
		int Tvertex;//List index of the terminate vertex in m_Svertex_list DE.
		get_integer(pDelimeter,pdstream,Tvertex);
		MGIges504Edge edge(curveDE,Svertex_list,Svertex,Tvertex_list,Tvertex);
		m_edges[i]=edge;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD504::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int num_edges=m_edges.size()-1;
	put_integer(num_edges,gsec,plines);
	for(int i=1; i<=num_edges; i++){
		const MGIges504Edge& edgei=m_edges[i];
		put_DEpointer(edgei.m_curve_DE,gsec,plines);
		put_DEpointer(edgei.m_Svertex_list,gsec,plines);
		put_integer(edgei.m_Svertex,gsec,plines);
		put_DEpointer(edgei.m_Tvertex_list,gsec,plines);
		put_integer(edgei.m_Tvertex,gsec,plines);
	}
}

//Get the start vertex.
MGBVertex* MGIges504Edge::get_SVertex(MGIgesIfstream& ifs)const{
	return ifs.m_vertexListMap.get_BVertex(m_Svertex_list,m_Svertex);
}

//Get the end vertex.
MGBVertex* MGIges504Edge::get_TVertex(MGIgesIfstream& ifs)const{
	return ifs.m_vertexListMap.get_BVertex(m_Tvertex_list,m_Tvertex);
}
