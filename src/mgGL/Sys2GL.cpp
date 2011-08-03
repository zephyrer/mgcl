/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// Sys2GL.cpp: mgSys2GL クラスのインプリメンテーション
#include "MGCLStdAfx.h"
#include <iostream>
#include "mgGL/Sys2GL.h"

// mgSys2GL

mgSys2GL::mgSys2GL(
   	size_t function_code,	//Function code.
	const MGGel* gel1, const MGGel* gel2
):mgSysGL(function_code,gel1), m_gel2(gel2){;}

//Test if this mgSysGL includes gel(return true) or not.
bool mgSys2GL::includes(const MGGel* gel)const{
	if(gel==object_id())
		return true;
	if(m_gel2==gel)
		return true;
	return false;
}

//replace gel_old to gel_new.
//If gel_old is not included in this, do nothing.
void mgSys2GL::replace(
	const MGGel* gel_old,	//gel_old must be a MGCurve.
	const MGGel* gel_new	//gel_new must be a MGFSurface.
){
	if(gel_old==object_id()){
		set_object_id(const_cast<MGGel*>(gel_new));
	}else{
		if(m_gel2==gel_old){
			m_gel2=gel_new;
		}
	}
}

// Output virtual function.
//Output to stream file:メンバデータを標準出力に出力する。
std::ostream& mgSys2GL::out(std::ostream& ostrm) const{
	mgSysGL::out(ostrm);
	ostrm<<",m_gel2="<<m_gel2;
	return ostrm;
}
