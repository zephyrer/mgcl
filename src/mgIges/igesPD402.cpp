/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD402.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Group.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD402.h"
#include "mgiges/IgesGsec.h"
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD402 is the class for Iges parameter data type 402(Group associativity).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD402.
MGIgesPD402::MGIgesPD402(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(ASSOCIATIVITY_INSTANCE,DEpointer){
}

//Read in parameter data from string stream data.
void MGIgesPD402::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int n;
	get_integer(pDelimeter,pdstream,n);//Number of entries.
	for(int i=0; i<n; i++){
		int entry;
		get_DEpointer(pDelimeter,pdstream,entry);
		m_members.push_back(entry);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD402::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int n=m_members.size();
	put_integer(n,gsec,plines);
	for(int i=0;i<n;i++)
		put_DEpointer(m_members[i],gsec,plines);
}

//Convert de(type=402: Associatibity entries) to MGGroup.
//Accepted type numbers are only group associativity: form number 1,7,14, and 15
//Returned is a newed MGGroup.
MGGroup* MGIgesIfstream::convert_group(
	const MGIgesDirectoryEntry& de
)const{
	int formn=de.FormNumber();
	if(formn!=1 && formn!=7 && formn!=14 && formn!=15)
		return 0;//Form number 1,7,14,or 15 is accepted.

	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD402* pd402=static_cast<const MGIgesPD402*>(pd.get());
	const std::vector<int>& members=pd402->m_members;
	int n=members.size();
	std::auto_ptr<MGGroup> group(new MGGroup);
	for(int i=0; i<n; i++){
		MGGel* geli=convert_to_gel(members[i]);
		if(geli)
			group->push_back(geli);
	}
	if(group->size())
		return group.release();
	else
		return 0;
}
