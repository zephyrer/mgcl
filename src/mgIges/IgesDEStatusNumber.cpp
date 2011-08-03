/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesDEStatusNumber.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesDEStatusNumber.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
using namespace std;

//!	@brief MGIgesDEStatusNumber describes the Status Number of a directory entry section.

// Constructors.
//! Constructs an object of class MGIgesDEStatusNumber.
//Default constructor, includes all the defalut value of MGCL.
MGIgesDEStatusNumber::MGIgesDEStatusNumber()
:m_BlankStatus(0), m_SubordinateEntitySwitch(independent),
m_EntityUseFlag(0), m_Hierarchy(0){
}

MGIgesDEStatusNumber::MGIgesDEStatusNumber(
	short BlankStatus,//0:Visible, 1:Blanked;
	SESwitch SubordinateEntitySwitch,
		//0:Independent, 1:Physically Dependent, 2: Logically Dependent, 3:Both 1 and 2.
	short EntityUseFlag,
		//0:Geometry, 1:Annotation, 2:Definition, 3:Other,
		//4:Logical/Positional, 5:2D Parametric, 6:Construction Geometry;
	short Hierarchy//0:Global top down, 1:Global defer, 2:Use hierarchy property;
):m_BlankStatus(BlankStatus),m_SubordinateEntitySwitch(SubordinateEntitySwitch),
m_EntityUseFlag(EntityUseFlag),m_Hierarchy(Hierarchy){
}

void MGIgesDEStatusNumber::read_in(const std::string& status){
	istringstream(status.substr(0,2))>>skipws>>m_BlankStatus;
	istringstream(status.substr(2,2))>>skipws>>m_SubordinateEntitySwitch;
	istringstream(status.substr(4,2))>>skipws>>m_EntityUseFlag;
	istringstream(status.substr(6,2))>>skipws>>m_Hierarchy;
}
