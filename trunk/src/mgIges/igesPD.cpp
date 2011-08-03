/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgIges/IgesPD.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

// Constructors.

//! Constructs an object of class MGIgesPDSec.
MGIgesPD::MGIgesPD()
:m_type_number(0),m_DEpointer(0){;}

//! Constructs an object of class MGIgesPDSec.
MGIgesPD::MGIgesPD(int type_number, MGIgesDirectoryEntry* DEpointer)
:m_type_number(type_number),m_DEpointer(DEpointer){;}

//Destructor;
MGIgesPD::~MGIgesPD(){;}
