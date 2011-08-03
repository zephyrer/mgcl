/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/DrawFunc.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/SPointSeq.h"
#include "mg/Group.h"
#include "mg/GelPositions.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/GLDrawFunc.h"
#include "mgGL/SysGLList.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Attach to MGCommandDrawer.
//Function's return value is the number of drawer's attached after attached.
size_t MGOpenGLView::attach_drawer(MGCommandDrawer* drawer){
	m_command_drawers.push_back(drawer);
	return m_command_drawers.size();
}

//Detach and Return the detached MGCommandDrawer.
void MGOpenGLView::detach_drawer(MGCommandDrawer* drawer){
	m_command_drawers.remove(drawer);
}

//////display member function.
void MGCurve::display_arrows()const{mgGDL::MGDrawArrow(*this,2);}
void MGCurve::display_break_points()const{
	const MGKnotVector& t=knot_vector();
	double ts=param_s(), te=param_e();
	size_t k=t.order(), n=t.bdim();
	for(size_t i=k-1; i<=n; i++){
		double tau=t[i];
		if(i>=k && tau==t[i-1]) continue;
		if(tau<ts || tau>te) continue;
		MGVector P=eval(tau);
		mgGDL::MGDrawPoint(P[0],P[1],P[2]);
	}
}
void MGLBRep::display_control_polygon()const{
	const MGBPointSeq& bp=line_bcoef();
	mgGDL::MGDrawPointSeq(bp);
}
void MGRLBRep::display_control_polygon()const{
	MGBPointSeq bp = non_homogeneous_bcoef();
	mgGDL::MGDrawPointSeq(bp);
}

void MGCurve::display_curvatures(
	double	scale,	//scaling of the graph.
	int		density,//densitiy of the graph.
	bool	use_radius//true:radius display, false:curvature display.
)const{
	mgGDL::MGDrawCurvaGraph(*this,scale,density,use_radius);
}

//////display member function.
void MGSurface::display_arrows()const{
	mgGDL::MGDrawArrow(*this,4,4);
}

//Display control polygons using mgGDL::MGDrawPointSeq(sp)
void MGSBRep::display_control_polygon()const{
	const MGSPointSeq& sp=surface_bcoef();
	mgGDL::MGDrawPointSeq(sp);
}

void MGRSBRep::display_control_polygon()const{
	MGSPointSeq sp=non_homogeneous_bcoef();
	mgGDL::MGDrawPointSeq(sp);
}

//////display member function.
void MGFace::display_arrows()const{mgGDL::MGDrawArrow(*this,4,4);}
void MGFace::display_control_polygon()const{
	surface()->display_control_polygon();
}

//////display member function.
void MGShell::display_arrows()const{
	size_t n=number_of_faces();
	for(size_t i=0; i<n; i++)
		face(i)->display_arrows();
}
void MGShell::display_control_polygon()const{
	size_t n=number_of_faces();
	for(size_t i=0; i<n; i++)
		face(i)->display_control_polygon();
}
