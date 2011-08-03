/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesParamDELine.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesParamLine.h"

//!	@brief MGIgesParamLine describes a line of Parameter Data of an IGES file.
//This is used to output one line os Parameter Data Section.

//! Constructs an object of class MGIgesParamLine.
MGIgesParamLine::MGIgesParamLine()
:m_DE_back_pointer(0){;}

//Construct inputting only DEpoiter.
MGIgesParamLine::MGIgesParamLine(int DEpointer)
:m_DE_back_pointer(DEpointer){;}

//one_line.size() must be <=64.
MGIgesParamLine::MGIgesParamLine(std::auto_ptr<std::string> one_line, int DEpointer)
:m_paramLine(one_line),m_DE_back_pointer(DEpointer){;};

