/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

//! @file
//!	@brief  Declaration for class MGIgesVertexListMap.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgIges/igesVertexListMap.h"
#include "mgIges/Igesifstream.h"
using namespace std;

//!	@brief MGIgesVertexListMap is the class to store MGBVertex*(newed objects) that
//are generated for MGIges502 vertices list.

//Obtain MGBVertex* of the the DE of the VERTEX List Entry(MGIgesPD502) and 
//the index vertex.
MGBVertex* MGIgesVertexListMap::get_BVertex(
	int DEid,	//Directory entry id of VERTEX list entry
	int vertex	  //List index of the VERTEX List Entry DE vertex_list
){
	const MGIgesDirectoryEntry* de=m_ifstream->directoryEntry(DEid);
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD502& pd502=*(static_cast<const MGIgesPD502*>(pd.get()));//reference to MGIgesPD502.
	int vID=m_VerticesVector.size();
	pair<PD502toVertexMap::iterator, bool> insertR=
		m_PD502toVertexMap.insert(make_pair(&pd502,vID));
	if(insertR.second){
		//If the pair(&pd502,vID) is inserted, make MGPvector<MGBVertex> for this pd502.
		size_t nvert=pd502.m_vertices.size();
		m_VerticesVector.push_back(std::vector<MGBVertex*>(nvert,0));//insert dummy array.
	}else
		//If not.
		vID=insertR.first->second;//Already inserted &pd502's array index.

	assert(vID<int(m_VerticesVector.size()));
	std::vector<MGBVertex*>& vertices=m_VerticesVector[vID];
	assert(vertex<int(vertices.size()));
	MGBVertex* vertexP=vertices[vertex];
	if(!vertexP){
		//This is the 1st refernece to pd502's vertex. Generate MGBVertex*.
		vertexP=new MGBVertex(pd502[vertex]);
		vertices[vertex]=vertexP;
	}
	return vertexP;//Return the found binder vertex.
}
