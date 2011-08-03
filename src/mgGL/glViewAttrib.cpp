/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <bitset>
#include "mgGL/OpenGLView.h"
#include "mgGL/glViewAttrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//////////////////////////////METHOD//////////////////

MGglViewAttrib::MGglViewAttrib(
	const MGOpenGLView& glview
):m_fovy(glview.m_fovy),
m_perspective(glview.m_perspective), m_near(glview.m_near),
m_far(glview.m_far),m_eyeP(glview.m_eyeP), m_up_vector(glview.m_up_vector),
m_center(glview.m_center),m_diameter(glview.m_diameter),m_scale(glview.m_scale),
m_cx(glview.m_cx), m_cy(glview.m_cy),m_cplane(glview.m_cplane){
	glview.make_RC_current();// cout<<"saving cplane:"<<m_cplane<<endl;
	glGetDoublev(GL_MODELVIEW_MATRIX,m_currentMat);
}

//Copy the informations of glview into this.
void MGglViewAttrib::copy(const MGOpenGLView& glview){
	m_fovy=glview.m_fovy;
	m_perspective=glview.m_perspective;
	m_near=glview.m_near; m_far=glview.m_far;
	m_eyeP=glview.m_eyeP;
	m_up_vector=glview.m_up_vector;
	m_center=glview.m_center; m_diameter=glview.m_diameter;
	m_scale=glview.m_scale;
	m_cx=glview.m_cx; m_cy=glview.m_cy;
	m_cplane=glview.m_cplane;

	glview.make_RC_current();
	glGetDoublev(GL_MODELVIEW_MATRIX,m_currentMat);
}

//set the viewing environment info.
void MGglViewAttrib::set(
	double neard, double fard,	//near and far data
	double diameter, const MGPosition& center,	//Model's size(diameter) and the center
	double scale, double cx, double cy,	//Current viewing's scale, and panning data(cx, cy).
	const MGPosition& eye, const MGVector& up_vector,//viewing eye position and view up vector
								//These parameters are used fro gluLookAt. see the explanation.
	const double currentMat[16]	//OpenGL's ModelView current matrix.
){
	m_near=neard; m_far=fard;
	m_diameter=diameter; m_center=center;
	m_scale=scale; m_cx=cx; m_cy=cy;
	m_eyeP=eye; m_up_vector=up_vector;
	for(size_t i=0; i<16; i++) m_currentMat[i]=currentMat[i];
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGglViewAttrib& atr){
	std::bitset<32> boolData;
	boolData.set(0,atr.m_perspective);

	buf << boolData.to_ulong();
	buf<<atr.m_fovy;
	buf<<atr.m_near;
	buf<<atr.m_far;
	atr.m_eyeP.dump(buf);
	atr.m_up_vector.dump(buf);
	atr.m_center.dump(buf);
	buf<<atr.m_diameter;
	buf<<atr.m_scale;
	buf<<atr.m_cx;
	buf<<atr.m_cy;
	buf<<atr.m_cplane;
	for(size_t i=0; i<16; i++) buf<<atr.m_currentMat[i];
	return buf;
}

MGIfstream& operator>> (MGIfstream& buf, MGglViewAttrib& atr){
	unsigned long lbit;
	buf >> lbit;
	std::bitset<32> boolData(lbit);
	atr.m_perspective=boolData[0];

	buf>>atr.m_fovy;
	buf>>atr.m_near;
	buf>>atr.m_far;
	atr.m_eyeP.restore(buf);
	atr.m_up_vector.restore(buf);
	atr.m_center.restore(buf);
	buf>>atr.m_diameter;
	buf>>atr.m_scale;
	buf>>atr.m_cx;
	buf>>atr.m_cy;
	buf>>atr.m_cplane;
	for(size_t i=0; i<16; i++) buf>>atr.m_currentMat[i];
	return buf;
}
	
//Debug Function.
std::ostream& operator<< (std::ostream& out, const MGglViewAttrib& atr){
	out<<std::endl<<"MGglViewAttrib="<<&atr<<",m_perspective="<<atr.m_perspective;
	out<<",m_fovy="<<atr.m_fovy<<",m_near="<<atr.m_near<<",m_far="<<atr.m_far;
	out<<",m_eyeP="<<atr.m_eyeP<<",m_up_vector="<<atr.m_up_vector<<",m_center="<<atr.m_center;
	out<<",m_diameter="<<atr.m_diameter<<",m_scale="<<atr.m_scale;
	out<<",(m_cx,m_cy)=("<<atr.m_cx<<","<<atr.m_cy<<")";
	out<<",m_cplane="<<atr.m_cplane<<std::endl;
	out<<",m_currentMat=["<<atr.m_currentMat[0];
	for(size_t i=1;i<16; i++) out<<","<<atr.m_currentMat[i];
	out<<"]";
	return out;
}
