/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD186.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Shell.h"
#include "mg/Group.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD186.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;

//!	@brief MGIgesPD186 is the class for Iges parameter data type 186
//(MSBO:Manifold Solid B-Rep Object Entity).

// Constructors.

//Read in parameter data from string stream data.
void MGIgesPD186::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_shell_DE);
	int orientation;
	get_integer(pDelimeter,pdstream,orientation);
	m_orientation=orientation?true:false;

	int num_voids;
	get_integer(pDelimeter,pdstream,num_voids);
	m_void_shells.resize(num_voids);
	m_orientations.resize(num_voids);
	for(int i=0; i<num_voids; i++){
		int shellDE;
		get_DEpointer(pDelimeter,pdstream,shellDE);
		m_void_shells[i]=shellDE;
		int orientation;
		get_integer(pDelimeter,pdstream,orientation);
		m_orientations[i]=orientation ? true:false;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD186::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_shell_DE,gsec,plines);
	int orientation=m_orientation ? 1:0;
	put_integer(orientation,gsec,plines);

	int num_voids=m_void_shells.size();
	put_integer(num_voids,gsec,plines);
	for(int i=0; i<num_voids; i++){
		put_DEpointer(m_void_shells[i],gsec,plines);
		int orientation=m_orientations[i] ? 1:0;
		put_integer(orientation,gsec,plines);
	}
}

//Convert de(type=186: Manifold Solid B-Rep Object) to MGGell.
//Returned is a newed MGGel.
//When MSBO has one or more void shells, MSBO is represented as a MGGroup.
//When MSBO does not have a void, MSBO is represented as a MGShell.
MGGel* MGIgesIfstream::convert_MSBO(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD186* pd186=static_cast<const MGIgesPD186*>(pd.get());

	MGIgesDirectoryEntry& deshell=*(m_DirectoryEntries[pd186->m_shell_DE]);
	MGShell* s0=convert_shell(deshell);
	if(!pd186->m_orientation)
		s0->negate();

	int nvoid=pd186->m_void_shells.size();
	if(nvoid==0)
		return s0;

	MGGroup* grp=new MGGroup(s0);
	for(int i=1; i<nvoid; i++){
		MGIgesDirectoryEntry& dei=*(m_DirectoryEntries[pd186->m_void_shells[i]]);
		MGShell* si=convert_shell(dei);
		if(!pd186->m_orientations[i])
			si->negate();
		//std::cout<<*f<<std::endl;;///////*********
		grp->push_back(si);
	}
	return grp;
}
