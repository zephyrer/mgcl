/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD504_H__)
#define __MGIGESPD504_H__

#include <vector>
#include "mg/Position.h"
#include "mgiges/IgesPD.h"
class MGBVertex;
class MGIgesIfstream;

///@cond
class MGIges504Edge{
public:

MGIges504Edge(
):m_curve_DE(0),m_Svertex_list(0),m_Svertex(0),m_Tvertex_list(0),m_Tvertex(0){;};

MGIges504Edge(
	int curve_DE,	///<DE pointer of the model space curve of MGIgesPD504.
	int Svertex_list,
		///<pointer to the DE of the VERTEX List Entry(MGPD502) for the start vertex.
	int Svertex,  ///<List index of the start vertex in m_Svertex_list DE.
	int Tvertex_list,
		///<pointer to the DE of the VERTEX List Entry(MGPD502) for the terminate vertex.
	int Tvertex  ///<List index of the terminate vertex in m_Svertex_list DE.
):m_curve_DE(curve_DE),m_Svertex_list(Svertex_list),m_Svertex(Svertex)
,m_Tvertex_list(Tvertex_list),m_Tvertex(Tvertex){;};

///Get the start vertex.
MGBVertex* get_SVertex(MGIgesIfstream& ifs)const;

///Get the end vertex.
MGBVertex* get_TVertex(MGIgesIfstream& ifs)const;

public:
	int m_curve_DE;	///<DE pointer of the model space curve of MGIgesPD504.
	int m_Svertex_list;
		///<pointer to the DE of the VERTEX List Entry(MGPD502) for the start vertex.
	int m_Svertex;  ///<List index of the start vertex in m_Svertex_list DE.
	int m_Tvertex_list;
		///<pointer to the DE of the VERTEX List Entry(MGPD502) for the terminate vertex.
	int m_Tvertex;  ///<List index of the terminate vertex in m_Svertex_list DE.
};
///@endcond

///MGIgesPD504 is the class for the Iges parameter data type 504(EDGE list) form 1.
///Type 504 data is a constituent of type 508(LOOP).
class MGIgesPD504: public MGIgesPD{
public:
	// Constructors.

/// Constructs an object of class MGIgesPD504.
MGIgesPD504(MGIgesDirectoryEntry* DEpointer=0);

///Destructor;
~MGIgesPD504(){;};

MGIges504Edge& operator[](int i){return m_edges[i];};
const MGIges504Edge& operator[](int i)const{return m_edges[i];};

///append an edge.
void push_back(const MGIges504Edge& edge){m_edges.push_back(edge);};

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

//Member data. These are set as public.
	///Vector of MGIges504Edge's.
	std::vector<MGIges504Edge> m_edges;///<m_edges[0] is dummy. This is because list index
			///<of PD508 starts from 1.
};

#endif /// __MGIGESPD504_H__
