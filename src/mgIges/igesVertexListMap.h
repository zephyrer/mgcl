/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESVertexListMap_H__)
#define __MGIGESVertexListMap_H__

#include <map>
#include <vector>
#include "mg/Pvector.h"
#include "topo/BVertex.h"
#include "mgIges/igesPD502.h"
class MGIgesIfstream;

///MGIgesVertexListMap is the class to store MGBVertex*(newed objects) that
///are generated for MGIges502 Vertices list.
class MGIgesVertexListMap{
public:

	typedef std::map<const MGIgesPD502*,int> PD502toVertexMap;
		//MGIgesPD502 is a vactor of vertices. So the 2nd int(say i) is
		//a index of m_VerticesVector. Thus m_VerticesVector[i] is
		//the vector of vertices of MGIgesPD502.

	/// Constructors.
	/// Constructs an object of class MGIgesVertexListMap.
	MGIgesVertexListMap(MGIgesIfstream* ifs=0):m_ifstream(ifs){;};

	///Obtain MGBVertex* of the the DE of the Vertices List Entry(MGIgesPD504) and 
	///the indgex vertex.
	MGBVertex* get_BVertex(
		int DEid,	///<Directory entry id of VERTEX list entry
		int vertex	  ///<List index of the VERTEX List Entry DE vertices_list
	);

	void set_ifstream(MGIgesIfstream* ifs){m_ifstream=ifs;};

//Member data.
private:

	MGIgesIfstream* m_ifstream;///<Iges Input stream pointer.

	///Map to get the array index of m_VerticesVector from the key MGIges502 pointer.
	///1st int is MGIges502 pointer, and 2nd is the array index of m_VerticesVector.
	PD502toVertexMap m_PD502toVertexMap;

	///m_VerticesVector[i] is a vector of MGBVertex*(newed objects) that are generated for
	///a type 502 DE, where i=m_PD502toVertexMap[MGIgesPD502*]. Let pd502 is a MGIgesPD502*,
	///and k is vetex list index of the type502 vertiex list,
	///then MGBVertex* vertex=m_VerticesVector[i][k], where i=m_PD502toVertexMap[pd502].
	std::vector< std::vector<MGBVertex*> > m_VerticesVector;
};

#endif /// __MGIGESVertexListMap_H__
