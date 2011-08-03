/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/Point.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGGeometry
// Implementation of MGGeometry.

// コンストラクタ
// 初期化なしでオブジェクトを作成する。
MGGeometry::MGGeometry():m_box(0){;}

//Copy constructor.
MGGeometry::MGGeometry(const MGGeometry& geo2):MGObject(geo2){
	if(geo2.m_box) m_box=new MGBox(*(geo2.m_box));
	else m_box=0;
}

//デストラクタ
MGGeometry::~MGGeometry(){
	if(m_box) delete m_box;
}

//Assignment.
MGGeometry& MGGeometry::set_geometry(const MGGeometry& geo2){
	MGObject::operator=(geo2);
	if(m_box){
		if(geo2.m_box)
			(*m_box)=*(geo2.m_box);
		else{
			delete m_box;
			m_box=0;
		}
	}else{
		if(geo2.m_box)
			m_box=new MGBox(*(geo2.m_box));
	}
	return *this;
}

//Return minimum box that includes whole of the geometry.
const MGBox& MGGeometry::box() const{
	if(!m_box) m_box=compute_box();
	return *m_box;
}

//Compute direction unit vector of the geometry.
MGUnit_vector MGGeometry::direction(const MGPosition& param) const{
	return mgZ_UVEC;
}

//Error allowed in the parameter space of the geometry.
double MGGeometry::parameter_error() const{
	MGBox prange=parameter_range();
	size_t n=prange.sdim();
	double error=0., ework;
	for(size_t i=0; i<n; i++){
		ework=prange(i).relative_error();
		error+=ework*ework;
	}
	return sqrt(error);
}
