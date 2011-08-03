/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

//! @file
//!	@brief  Declaration for class MGIges504EdgeListMap.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgIges/IgesVertexListMap.h"
#include "mgIges/iges504EdgeListMap.h"
#include "mgIges/Igesifstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace std;

//!	@brief MGIges504EdgeListMap is the class to store MGEdge*(newed objects) that
//are generated for MGIges504 EDGE list.

//Obtain MGEdge* of the the DE of the EDGE List Entry(MGIgesPD504) and 
//the indgex edge.
MGEdge* MGIges504EdgeListMap::get_egde(
	int DEid,	//DE id of type 504
	int edge	//List index of the EDGE List Entryt DE edge_list
){
	const MGIgesDirectoryEntry* de=m_ifstream->directoryEntry(DEid);
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD504& pd504=*(static_cast<const MGIgesPD504*>(pd.get()));//reference to MGIgesPD504.
	int vID=m_EdgesVector.size();
	pair<PD504toEdgesMap::iterator, bool> insertR=
		m_PD504PtoVecID.insert(make_pair(&pd504,vID));
	if(insertR.second){
		//If the pair(&pd504,vID) is inserted, make MGPvector<MGEdge> for this pd504.
		size_t nedges=pd504.m_edges.size();
		m_EdgesVector.push_back(std::vector<MGEdge*>(nedges,0));//insert dummy array.
	}else//If not.
		vID=insertR.first->second;//Already inserted &pd504's array index.

	assert(vID<int(m_EdgesVector.size()));
	std::vector<MGEdge*>& edges=m_EdgesVector[vID];
	assert(edge<int(edges.size()));
	MGEdge* edgeP=edges[edge];
	if(!edgeP){
		//This is the 1st refernece to pd504's edge. Generate MGEdge*.
		const MGIges504Edge& edge504=pd504[edge];
		MGCurve* crv=static_cast<MGCurve*>(m_ifstream->convert_to_gel(edge504.m_curve_DE));
		if(crv){
			edgeP=new MGEdge(crv);
			MGBVertex* bvS=edge504.get_SVertex(*m_ifstream);
			MGBVertex* bvT=edge504.get_TVertex(*m_ifstream);
			edgeP->set_i_th_binder(0,*bvS);
			edgeP->set_i_th_binder(1,*bvT);
			edges[edge]=edgeP;
		}
	}
	return edgeP;//Return the found binder edge.
}
