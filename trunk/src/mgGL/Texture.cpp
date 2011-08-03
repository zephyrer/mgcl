/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Texture.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
	
//Invoke appropriate OpenGL fucntion to this attribute.
void MGTexture::exec()const{return;}

//Debug Function.
std::ostream& MGTexture::out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"MGTexture="<<this;
	return ostrm;
}

//Read all member data.
MGOfstream& operator<< (MGOfstream& buf, const MGTexture& txtr){
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGTexture& txtr){
	return buf;
}
