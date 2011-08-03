/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD508_H__)
#define __MGIGESPD508_H__

#include <vector>
#include "mg/Position.h"
#include "mgiges/IgesPD.h"
class MGEdge;
class MGBVertex;
class MGIgesIfstream;

///@cond
class MGIges508Edge{
public:

MGIges508Edge(
):m_type(0),m_edge_list(0),m_edge(0),m_orientation(0)
,m_isoparameterics(1),m_pcurves(1){;};

MGIges508Edge(
short type,	///<type of the edge, =0:edge, =1:vertex.
int edge_list,
	///<pointer to the DE of the VERTEX List Entry(MGPD502) for the start vertex.
int edge,  ///<List index of the edge into edge_list.
short orientation///<Orientation of the edge.
		///<=0: opposite, =1: the directions agree.
):m_type(type),m_edge_list(edge_list),m_edge(edge),m_orientation(orientation){;};

///Get the edge pointer(newed object).
MGEdge* get_edge(const MGIgesIfstream& ifs)const;

///Get the edge pointer(newed object).
MGBVertex* get_BVertex(const MGIgesIfstream& ifs)const;

bool is_vertex()const{return m_type!=0;};
bool direction_is_opposite()const{return m_orientation==0;};
int number_of_pcurves()const{return m_pcurves.size()-1;};

///Read in parameter data from string stream data.
void read_in(
	char pDelimeter,
	std::istringstream& pdstream
);

///Write out this MGIges508Edge as MGIgesParamLine's(into plines).
///Except for string data, one integer or double data is output
///into one MGIgesParamLine, not striding over more than one line.
///Only when string data is output(to Holleris string), the data
///may stride over more than one lines.
///plines[i] for 0<=i<plines.size() are valid.
void write_out_into_string(
	const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines ///<output plines.
)const;

public:
	short m_type;	///<type of the edge, =0:edge, =1:vertex.
	int m_edge_list;
		///<pointer to the DE of the EDGE(vertex) List Entry(MGIgesPD504/MGIgesPD502)
		///<for the edge.

	int m_edge;  ///<List index of the EDGE(vertex) List Entryt DE m_edge_list.
	short m_orientation;///<Orientation of the edge.
		///<=0: opposite, =1: the directions agree, to the model space curve.
	
	std::vector<bool> m_isoparameterics;///<vector of the bools that indicate
		///<whether m_pcurves[i] is isoparametric(m_isoparameterics[i]=true) or not.

	std::vector<int> m_pcurves;  ///<vector of the parameter space curve DE pointer.
	///<******* m_pcurves[0] and m_isoparameterics[0] are dummy. This is because igesPD510
	///<list index start from 1.
	///<m_pcurves.size()==m_isoparameterics.size()==number of edges +1.
};
///@endcond

///MGIgesPD508 is the class for Iges parameter data type 508(LOOP).
class MGIgesPD508: public MGIgesPD{
public:
	// Constructors.

	/// Constructs an object of class MGIgesPD508.
	MGIgesPD508(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD508(){;};

	MGIges508Edge& edge(size_t i){return *(m_edges[i]);};
	const MGIges508Edge& edge(size_t i)const{return *(m_edges[i]);};

	///Get the number of edges.
	size_t number_of_edges()const{return m_edges.size()-1;};

	///append an edge, which must be a newed one and whose ownership will be
	///transfered to this.
	void push_back(MGIges508Edge* edge){m_edges.push_back(edge);};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	///Write out this PD as MGIgesParamLine's(into plines).
	///Except for string data, one integer or double data is output
	///into one MGIgesParamLine, not striding over more than one line.
	///Only when string data is output(to Holleris string), the data
	///may stride over more than one lines.
	///plines[i] for 0<=i<plines.size() are valid.
	void write_out_into_string(
		const MGIgesGSec& gsec,	///Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///output plines.
	)const;

private:

//Member data. These are set as public.
	///Vertices of 3D coordinates.
	MGPvector<MGIges508Edge> m_edges;///m_edges[0] is dummy.
			///This is because list index of MGIges508Edge's m_edge starts from 1.
};

#endif // __MGIGESPD508_H__
