/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD502.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD502.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD502 is the class for Iges parameter data type 502(VERTEX list).

// Constructors.

//! Constructs an object of class MGIgesPD502.
MGIgesPD502::MGIgesPD502(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(VERTEX,DEpointer),m_vertices(1,MGPosition()){
}

//append one vertex data.
void MGIgesPD502::push_back(
	const MGPosition& vertex
){
	m_vertices.push_back(vertex);
}

//Read in parameter data from string stream data.
void MGIgesPD502::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_vertices;
	get_integer(pDelimeter,pdstream,num_vertices);
	m_vertices.resize(num_vertices+1);
	for(int i=1; i<=num_vertices; i++){
		MGPosition P(3);
		for(int j=0; j<3; j++)
			get_real(pDelimeter,pdstream,P(j));
		m_vertices[i]=P;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD502::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int num_vertices=m_vertices.size()-1;
	put_integer(num_vertices,gsec,plines);
	for(int i=1; i<=num_vertices; i++){
		const MGPosition& P=m_vertices[i];
		put_real(P[0],gsec,plines);
		put_real(P[1],gsec,plines);
		put_real(P[2],gsec,plines);
	}
}
