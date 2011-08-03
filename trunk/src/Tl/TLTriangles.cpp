/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "Tl/TLPoints.h"
#include "Tl/TLTriangles.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Return the i-th mgTLTriangle.
const mgTLTriangle& mgTLTriangles::operator[](size_t i)const{
	return *(m_triangles[i]);
}
mgTLTriangle& mgTLTriangles::operator[](size_t i){
	return *(m_triangles[i]);
}

std::ostream& operator<< (std::ostream& out, const mgTLTriangles& tlTriangles){
	out<<"TLTriangles="<<(&tlTriangles)<<" num of triangles="
		<<tlTriangles.m_triangles.size()<<std::endl;
	mgTLTriangles::const_triIterator iter = tlTriangles.begin();
	for(; iter != tlTriangles.end(); iter++){out<<**iter;}
	out<<std::endl;
	return out;
}

void mgTLTriangles::print(std::ostream& out, const mgTLData& tld)const{
	out<<"TLTriangles="<<this<<" num of triangles="
		<<this->m_triangles.size()<<std::endl;
	const_triIterator iter=this->begin(), iend=this->end();
	for(; iter != iend; iter++){(**iter).print(out,tld);}
	out<<std::endl;
}
