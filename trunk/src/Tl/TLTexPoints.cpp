/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "Tl/TLTexPoints.h"

//mgTLTexPoints holds the vector of texture coordinates of mgTLPoints.
//the surface parameter (u,v).
using namespace std;

////////// constructor /////////////
mgTLTexPoints::mgTLTexPoints(
	size_t n			//number of points.
	//double image_width,	//image width in world coordinates.
	//double image_height	//image height in world coordinates.
):m_tex_positions(n),m_refs(n){
//m_image_width(image_width),m_image_height(image_height){
	for(size_t i=0; i<n; i++){
		m_refs[i]=0;
	}
}

///////// member function ////////////

void mgTLTexPoints::set_texcoord(size_t i, double s, double t){
	MGPosition& st=m_tex_positions[i];
	int nref=m_refs[i];
	if(nref<0) return;
	//if(nref>0){
	//	st*=double(nref);
	//	st(0)+=s; st(1)+=t;
	//	st/=double(nref+1);
	//}else{
		st=MGPosition(s,t);
	//}
	m_refs[i]+=1;
}
void mgTLTexPoints::set_texcoord(size_t i, const MGPosition &texc){
	set_texcoord(i,texc[0],texc[1]);
}
void mgTLTexPoints::set_texcoord_as_not_modify(size_t i, double s, double t){
	m_tex_positions[i]=MGPosition(s,t);
	m_refs[i]=-1;
}
void mgTLTexPoints::set_texcoord_as_not_modify(size_t i, const MGPosition &texc){
	set_texcoord_as_not_modify(i,texc[0],texc[1]);
}

ostream& operator<< (ostream& out, const mgTLTexPoints& Points){
	size_t n=Points.size();
	out<<"TLTexPoints="<<(&Points)<<" num="<<n<<endl;
	for(size_t i=0; i<n; i++){
		out<<i<<":(s,t)="<<Points.m_tex_positions[i]<<", ref="<<Points.m_refs[i]<<endl;
	}
	return out;
}
