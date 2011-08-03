/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCellMap_HH_
#define _MGCellMap_HH_
#include <map>
#include "mg/MGCL.h"

///@cond

class MGCellNB;

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::map<const MGCellNB*, MGCellNB*>;
#pragma warning( pop )
#endif

class MGCLASS MGCellMap{
public:
	typedef std::map<const MGCellNB*, MGCellNB*> map;

	MGCellMap(){;};
	virtual ~MGCellMap(){;};

	std::map<const MGCellNB*, MGCellNB*> m_map;

};

///@endcond

#endif
