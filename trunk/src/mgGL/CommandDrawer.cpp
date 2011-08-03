/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <functional>
#include "mg/DrawFunc.h"
#include "mgGL/GLDrawFunc.h"
#include "mgGL/CommandDrawer.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// implementation for class MGCommandDrawer.

void MGCommandDrawer::draw_points(
	const MGColor& color,
	const std::vector<MGPosition>& ipos
)const{
	mgGDL::draw_points(color,MGDrawFunc::PColor2(),ipos);
}

void MGCommandDrawer::draw_lines(
	const MGColor& color,
	const std::vector<MGPosition>& ipos
)const{
	if(ipos.size()>=2){
		color.exec();//line color
		mgGDL::MGDrawPolyline(ipos);
	}
}

void MGCommandDrawer::draw_points_lines(
		const MGColor& pcolor,	//color of points(boundary)
		const MGColor& lcolor,	//color of line.
		const std::vector<MGPosition>& ipos
)const{
	draw_lines(lcolor,ipos);
	draw_points(pcolor,ipos);
}
