/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGES504EdgeListMap_H__)
#define __MGIGES504EdgeListMap_H__

#include <map>
#include <vector>
#include "mg/Pvector.h"
#include "topo/Edge.h"
#include "mgIges/igesPD504.h"
class MGIgesIfstream;

///	@brief MGIges504EdgeListMap is the class to store MGEdge*(newed objects) that
///are generated for MGIges504 EDGE list.
class MGIges504EdgeListMap{
public:

	///key is MGIgesPD504* and the data is the array index i of m_EdgesVector.
	///m_EdgesVector[i] is MGIgesPD504's data, MGPvector<MGEdge>.
	typedef std::map<const MGIgesPD504*,int> PD504toEdgesMap;

	/// Constructors.
	/// Constructs an object of class MGIges504EdgeListMap.
	MGIges504EdgeListMap(MGIgesIfstream* ifs=0):m_ifstream(ifs){;};

	///Obtain MGEdge* of the the DE of the EDGE List Entry(MGIgesPD504) and 
	///the indgex edge.
	MGEdge* get_egde(
		int DEid,	///<DE id of type 504
		int edge	 ///<List index of the EDGE List Entry DE edge_list
	);

	void set_ifstream(MGIgesIfstream* ifs){m_ifstream=ifs;};

////Member data.
private:

	MGIgesIfstream* m_ifstream;///Iges Input stream pointer.

	///Map to get the array index of m_EdgesVector from the key MGIges504 pointer.
	///1st is MGIges504 pointer, and 2nd is the array index of m_EdgesVector.
	PD504toEdgesMap m_PD504PtoVecID;

	///m_EdgesVector[i] is a vector of MGEdge*(newed objects) that are generated for
	///the DE j(or MGIgesPD504* pd504), where i=m_PD504PtoVecID[pd504].
	///Let k is the edge list index of MGIges504Edge DE j,
	///then MGEdge* edge=m_EdgesVector[i][k], where i=m_EdgeDEtoVecID[pd504].
	std::vector< std::vector<MGEdge*> > m_EdgesVector;
};

#endif // __MGIGES504EdgeListMap_H__
