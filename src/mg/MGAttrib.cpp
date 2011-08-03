/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Attrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mgGL/Context.h"
#include "mgGL/GLAttrib.h"
#include "mg/GelFactory.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGAttrib
// Implementation of MGAttrib.

//Virtual Destructor
MGAttrib::~MGAttrib(){;}

//Construct a null newed MGAttrib from the type id TID.
MGAttrib* MGNullAttrib(long TID){
	MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
	return static_cast<MGAttrib*>(reg->create_gel(TID));
}

AUTO_GEL_REGISTER(MGContext, MGCONTEXT_TID);

//Read all member data.
void MGAttrib::ReadMembers(MGIfstream& buf){;}
//Write all member data
void MGAttrib::WriteMembers(MGOfstream& buf)const{;}

