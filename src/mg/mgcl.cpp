/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mgGL/color.h"
#include "mg/DrawFunc.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of non-categolized global functions.
//

//Version no. definition. This is used in MGIfstream/MGOfstream.cpp.
MGEXTERN const char* _MGCL_VER="MGCL0800";
MGEXTERN const char* _MGCL_FILE="File DG Technologies,Inc";
MGColors MGColor::m_colors=MGColors();

const MGColor& MGDrawFunc::m_PColor1=MGColor::get_instance(MGColor::Black);//Boundary color of the point.
const MGColor& MGDrawFunc::m_PColor2=MGColor::get_instance(MGColor::White);//Inner square color of the point.
const float MGDrawFunc::m_point_size=5.;	//Point rectangle width and height of below.

//Get the MGCL_Version number.
const char* MGCL_Version(){
	return _MGCL_VER;
}
//Get the MGCL File validity.
const char* MGCL_File_validity(){
	return _MGCL_FILE;
}

namespace MGCL{

void start_up(bool need_to_GdiStartUp){
	if(need_to_GdiStartUp){
		//-> GDI+ startup
		Gdiplus::GdiplusStartupInput gdiplusStartupInput;
		Gdiplus::GdiplusStartup( &m_gdiplusToken, &gdiplusStartupInput, NULL );
		m_gdiplus_initialized=true;
		//<- GDI+ startup
	}else{
		m_gdiplus_initialized=false;
	}
}

void shut_down(){
	if(m_gdiplus_initialized)
		Gdiplus::GdiplusShutdown( m_gdiplusToken );
}

}
