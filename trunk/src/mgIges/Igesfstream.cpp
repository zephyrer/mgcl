/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

//! @file
//!	@brief  Declaration for class MGIgesFstream.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgIges/Igesfstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//Initialize all the member data to the state of no_value_holding.
void MGIgesFstream::initialize(const char* filename){
	m_StartSection=std::string();
	m_GSection=MGIgesGSec(filename);
	m_nlineGSec=0;
	m_DirectoryEntries.clear();
	MGIgesDirectoryEntry* de=new MGIgesDirectoryEntry;
	m_DirectoryEntries.push_back(de);//Set the dummy record.
}

//Function's return value is the directory entry pointer pushed back.
int MGIgesFstream::push_back_DE(MGIgesDirectoryEntry* de){
	int deNum=m_DirectoryEntries.size();
	m_DirectoryEntries.push_back(de);
	return deNum;
}
