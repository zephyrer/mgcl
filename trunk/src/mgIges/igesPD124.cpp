/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD124.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD124.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD124 is the class for Iges parameter data type 124(Transformation matrix).

// Constructors.

//! Constructs an object of class MGIgesPD124.
MGIgesPD124::MGIgesPD124(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(TRANSFORMATION_MATRIX,DEpointer){
	for(int i=0; i<12; i++)
		m_matrix[i]=0.;
}

MGIgesPD124::MGIgesPD124(const MGTransf& tr)
:MGIgesPD(TRANSFORMATION_MATRIX){
	m_matrix[0]=tr.ref(0,0);m_matrix[1]=tr.ref(1,0);m_matrix[2]=tr.ref(2,0); m_matrix[3]=tr.ref(3,0);
	m_matrix[4]=tr.ref(0,1);m_matrix[5]=tr.ref(1,1);m_matrix[6]=tr.ref(2,1); m_matrix[7]=tr.ref(3,1);
	m_matrix[8]=tr.ref(0,2);m_matrix[9]=tr.ref(1,2);m_matrix[10]=tr.ref(2,2);m_matrix[11]= tr.ref(3,2);
}

//convert this transformation to MGTransf.
void MGIgesPD124::convert_to_MGTransf(MGTransf& tr)const{
	tr.resize(3);
	tr(0,0)=m_matrix[0];tr(1,0)=m_matrix[1];tr(2,0)=m_matrix[2]; tr(3,0)=m_matrix[3];
	tr(0,1)=m_matrix[4];tr(1,1)=m_matrix[5];tr(2,1)=m_matrix[6]; tr(3,1)=m_matrix[7];
	tr(0,2)=m_matrix[8];tr(1,2)=m_matrix[9];tr(2,2)=m_matrix[10]; tr(3,2)=m_matrix[11];
}

//Read in parameter data from string stream data.
void MGIgesPD124::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	const int msize=12;//m_matrixÇÃóvëfêî
	for(int i=0; i<msize; i++)
		MGIges::get_real(pDelimeter,pdstream,m_matrix[i]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD124::write_out_into_string(
	const MGIgesGSec& gsec,	
	MGPvector<std::string>& plines 
)const{
	const int msize=12;//m_matrixÇÃóvëfêî
	for(int i=0;i<msize ;i++)
		MGIges::put_real(m_matrix[i],gsec,plines);
}

