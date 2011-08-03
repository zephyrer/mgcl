/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/AbstractGels.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGAbstractGels Class.

///////////////operator overloaded//////////////

//////////Member Function//////////

// Output virtual function.
std::ostream& operator<<(std::ostream& ostrm, const MGAbstractGels& agells){
	ostrm<<"MGAbstractGels::number of agells="<<agells.size()<<std::endl;
	MGAbstractGels::const_iterator i=agells.begin(), ie=agells.end();	
	for(size_t j=0; i!=ie; i++, j++){
		ostrm<<"agell-"<<j<<":"<<(*i)<<std::endl;
	}
	return ostrm;
}
