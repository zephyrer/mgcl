/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <bitset>
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/CSisect.h"
#include "mg/Straight.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/ConstructionPlane.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//construction plane's color definitions.
const static float lColorV[4]={0.6f, 0.6f, 0.6f, 1.f};//line color except (u,v) plus axis.

MGConstructionPlane::MGConstructionPlane(
	const MGPlane& plane,//construction plane.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
):m_disabled(false),m_bind_to_grid(false)
,m_plane(MGUnit_vector(plane.u_deriv()),MGUnit_vector(plane.v_deriv()),MGPosition(plane.root_point()))
,m_uspan(uspan), m_vspan(vspan), m_nspan(nspan)
,m_unum(uline_num), m_vnum(vline_num)
,m_uaxisColor(MGColor::get_instance(MGColor::Red))
,m_vaxisColor(MGColor::get_instance(MGColor::Green)){
	m_lineColor.set_color(lColorV);
}

MGConstructionPlane::MGConstructionPlane(
	double origin[3],	//Origin's coordinate value.
	double uaxis[3],	//A vector value of the horizontal direction.
	double vaxis[3],	//A vector value of the vertical direction.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
):m_disabled(false),m_bind_to_grid(false)
,m_uaxisColor(MGColor::get_instance(MGColor::Red))
,m_vaxisColor(MGColor::get_instance(MGColor::Green)){
	m_lineColor.set_color(lColorV);
	m_plane=MGPlane(
		MGUnit_vector(MGVector(3,uaxis)),MGUnit_vector(MGVector(3,vaxis)),MGPosition(3,origin));
}

//Construct a construction plane froma a box of 3D space.
MGConstructionPlane::MGConstructionPlane(
	const MGBox& box,
	int view_num	//Standard view number:
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
):m_disabled(false),m_bind_to_grid(false)
,m_uaxisColor(MGColor::get_instance(MGColor::Red))
,m_vaxisColor(MGColor::get_instance(MGColor::Green)){
	m_lineColor.set_color(lColorV);
	set_grid_data(box,view_num);
}

//Compute grid data and the plane from the plane and the grid span data.
void MGConstructionPlane::set_grid_data(
	const MGPlane& plane,//construction plane.
	double uspan,		//span length along u axis.
	double vspan,		//span length along v axis.
	int uline_num,		//number of lines along u axis.
	int vline_num,		//number of lines along v axis.
	double nspan		//span length along normal axis.
){
	m_plane=MGPlane(
		MGUnit_vector(plane.u_deriv()),MGUnit_vector(plane.v_deriv()),MGPosition(plane.root_point()));
	m_uspan=uspan;
	m_vspan=vspan;
	m_nspan=nspan;
	m_unum=uline_num;
	m_vnum=vline_num;
}

//Compute grid data and the plane from box and view.
void MGConstructionPlane::set_grid_data(
	const MGBox& box,
	int view_num	//Standard view number:
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
){
	double span;
	size_t lnum;
	MGPosition mid;

	size_t sdid;//maxmum area coordinate pair will be output.
				//0:(x,y), 1:(y,z), 2:(z,x)
	MGcplane_parameter(box,span,lnum,sdid,mid);
	if(view_num>=2 && view_num<=4)
		sdid=size_t(view_num-2);
	set_span(span);
	set_num(lnum);
	MGVector uderi(0.,0.,0.),vderi(0.,0.,0.);
	size_t sdid2=(sdid+2)%3, sdid1=(sdid+1)%3;
	uderi(sdid)=1.; vderi(sdid1)=1.;
	plane()=MGPlane(uderi,vderi,mid);
}

//Bind the input point uv(this construction plane's parameter value)
//to the nearest grid point of this construction plane.
void MGConstructionPlane::bind_to_grid(
	const MGPosition& uv, MGPosition& uvout
)const{
	uvout.resize(2);
	double x=double(int(fabs(uv[0])/m_uspan+.5))*m_uspan;
	if(uv[0]<0.) uvout(0)=-x; else uvout(0)=x;
	double y=double(int(fabs(uv[1])/m_vspan+.5))*m_vspan;
	if(uv[1]<0.) uvout(1)=-y; else uvout(1)=y;
}

//Convert cplane coordinates to the normal world coordinates.
//Function's return is the world coordinates.
MGPosition MGConstructionPlane::convert_to_world(const MGPosition& cplane_coord)const{
	MGPosition world(m_plane.root_point());
	world+=m_uspan*cplane_coord[0]*m_plane.u_deriv();
	world+=m_vspan*cplane_coord[1]*m_plane.v_deriv();
	world+=m_nspan*cplane_coord[2]*m_plane.normal();
	return world;
}

//Convert to cplane coordinates from the normal world coordinates.
//Function's return is the cplane coordinates.
MGPosition MGConstructionPlane::convert_from_world(const MGPosition& world_coord)const{
	MGPosition cplaneP(world_coord-m_plane.root_point());
	cplaneP(0)/=m_uspan;
	cplaneP(1)/=m_vspan;
	cplaneP(2)/=m_nspan;
	return cplaneP;
}

//Get line and axis colors
void MGConstructionPlane::get_color(
	MGColor& lineColor,		//Grid line color
	MGColor& uaxisColor,	//u axis
	MGColor& vaxisColor		//v axis
)const{
	lineColor=m_lineColor;
	uaxisColor=m_uaxisColor;
	vaxisColor=m_vaxisColor;
}

//locate a point on this plane, given straight line.
//the located point will be the intersection(or nearest bound) point
//of sl and the plane.
//Function's return value is the world coordinate located.
MGPosition MGConstructionPlane::locate(
	const MGStraight& sl,//input the ray straight line of the cursor.
	MGPosition& uv		//the plane's parameter value (u,v) will be output.
)const{
	MGCSisect is;
	sl.relation(plane(), is);
	if(is_bind_to_grid()){
		bind_to_grid(is.param_surface(),uv);
	}else uv=is.param_surface();
	return plane().eval(uv);
}

//Draw this plane using OpenGL.
void MGConstructionPlane::draw()const{
	if(disabled()||!valid())
		return;
	glPushAttrib(GL_LINE_BIT);//=========GL_ALL_ATTRIB_BITS
	glDisable(GL_LIGHTING);
	glDisable(GL_LINE_STIPPLE);

	const MGPosition& origin=m_plane.root_point();
	const double* origind=origin.data();
	MGVector udire=m_plane.u_deriv()*m_uspan;
	MGVector vdire=m_plane.v_deriv()*m_vspan;
	MGPosition u0(origin-udire*m_unum), u1(origin+udire*m_unum);
	MGPosition v0(origin-vdire*m_vnum), v1(origin+vdire*m_vnum);

	int i,j;
	MGPosition vS1(v0), vE1(v1), vS2(v0), vE2(v1);
	MGPosition uS1(u0), uE1(u1), uS2(u0), uE2(u1);

	glLineWidth(1.f);
	m_lineColor.exec();
	glBegin(GL_LINES);
		for(i=10; i<=m_unum; i+=10){
			for(int i1=1; i1<=9; i1++){
				vS1+=udire; vE1+=udire;
				glVertex3dv(vS1.data());glVertex3dv(vE1.data());
				vS2-=udire; vE2-=udire;
				glVertex3dv(vS2.data());glVertex3dv(vE2.data());
			}
			vS1+=udire; vE1+=udire;
			vS2-=udire; vE2-=udire;
		}
		for(j=10; j<=m_vnum; j+=10){
			for(int j1=1; j1<=9; j1++){
				uS1+=vdire; uE1+=vdire;
				glVertex3dv(uS1.data());glVertex3dv(uE1.data());
				uS2-=vdire; uE2-=vdire;
				glVertex3dv(uS2.data());glVertex3dv(uE2.data());
			}
			uS1+=vdire; uE1+=vdire;
			uS2-=vdire; uE2-=vdire;
		}
	glEnd();

	glLineWidth(2.f);
	////u minus axis and v minus axis.
	glBegin(GL_LINE_STRIP);
		glVertex3dv(u0.data());glVertex3dv(origind);glVertex3dv(v0.data());
	glEnd();

	vS1=vS2=v0; vE1=vE2=v1;
	uS1=uS2=u0; uE1=uE2=u1;
	glBegin(GL_LINES);
		MGVector udire10=udire*10.;
		for(i=10; i<=m_unum; i+=10){
			vS1+=udire10; vE1+=udire10;
			glVertex3dv(vS1.data());glVertex3dv(vE1.data());
			vS2-=udire10; vE2-=udire10;
			glVertex3dv(vS2.data());glVertex3dv(vE2.data());
		}
		MGVector vdire10=vdire*10.;
		for(j=10; j<=m_vnum; j+=10){
			uS1+=vdire10; uE1+=vdire10;
			glVertex3dv(uS1.data());glVertex3dv(uE1.data());
			uS2-=vdire10; uE2-=vdire10;
			glVertex3dv(uS2.data());glVertex3dv(uE2.data());
		}
	glEnd();

	glLineWidth(2.f);
	m_uaxisColor.exec();////u plus axis.
	glBegin(GL_LINES);
		glVertex3dv(origind);glVertex3dv(u1.data());
	glEnd();
	m_vaxisColor.exec();////v plus axis.
	glBegin(GL_LINES);
		glVertex3dv(origind);glVertex3dv(v1.data());
	glEnd();
	glPopAttrib();
}

//set the colors to default ones.
void MGConstructionPlane::set_default_color(){
	m_uaxisColor=MGColor::get_instance(MGColor::Red);
	m_vaxisColor=MGColor::get_instance(MGColor::Green);
	m_lineColor.set_color(lColorV);
}

//set grid line color.
void MGConstructionPlane::set_line_color(const MGColor& color){
	m_lineColor=color;
}

//set uaxis color.
void MGConstructionPlane::set_uaxis_color(const MGColor& color){
	m_uaxisColor=color;
}

//set vaxis color.
void MGConstructionPlane::set_vaxis_color(const MGColor& color){
	m_vaxisColor=color;
}

#define MESH_NUM 500.
//Compute cplane parameter from 3D box.
void MGcplane_parameter(
	const MGBox& box,
	double& span,	//span length will be output.
	size_t& lnum,	//number of lines along vertical and horizontal will be output.
	size_t& sdid,	//maxmum area coordinate pair will be output.
					//0:(x,y), 1:(y,z), 2:(z,x)
	MGPosition& mid)//rounded mid point will be output.
{
	double len[3];
	mid.resize(3);
	size_t i;
	for(i=0; i<3; i++){
		const MGInterval& intrvl=box[i];
		len[i]=intrvl.length().value();
		mid(i)=intrvl.mid_point();
	}
	double maxlen=len[0];
	if(maxlen<len[1]){
		if(len[1]<len[2]) maxlen=len[2];
		else maxlen=len[1];
	}else{
		if(maxlen<len[2]) maxlen=len[2];
	}
	span=maxlen*1.1/MESH_NUM;
	if(span>10000.) span=100000.;
	else if(span>1000.) span=10000.;
	else if(span>100.) span=1000.;
	else if(span>10.) span=100.;
	else if(span>1.) span=10.;
	else if(span>.1) span=1.;
	else if(span>.01) span=.1;
	else span=.01;
	lnum=size_t(MESH_NUM*.5);

	double area0=len[0]*len[1], area1=len[1]*len[2], area2=len[2]*len[0];
	sdid=0;
	if(area0<area1){
		if(area1<area2) sdid=2;
		else sdid=1;
	}else{
		if(area0<area2) sdid=2;
	}
//	double span10=span*10.;
//	double hspan10=span10*.5;
	double hspan=span*.5;
	for(i=0; i<3; i++)
		mid(i)=int((mid(i)+hspan)/span)*span;
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGConstructionPlane& cpl){
	std::bitset<32> boolData;
	boolData.set(0,cpl.m_disabled);
	boolData.set(1,cpl.m_bind_to_grid);

	buf << boolData.to_ulong();
	buf.WritePointer(&(cpl.m_plane));
	buf << cpl.m_vspan;
	buf << cpl.m_uspan;
	buf << cpl.m_vnum;
	buf << cpl.m_unum;
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGConstructionPlane& cpl){
	unsigned long boollong;
	buf >> boollong;
	std::bitset<32> boolData(boollong);
	cpl.m_disabled=boolData[0];
	cpl.m_bind_to_grid=boolData[1];

	MGGel* gel=	buf.ReadPointer();
	MGPlane* pl=static_cast<MGPlane*>(gel);
	cpl.m_plane=*pl; delete pl;
	buf >> cpl.m_vspan;
	buf >> cpl.m_uspan;
	buf >> cpl.m_vnum;
	buf >> cpl.m_unum;
	return buf;
}
	
//Debug Function.
std::ostream& operator<< (std::ostream& out, const MGConstructionPlane& pln){
	out<<std::endl<<"MGConstructionPlane="<<&pln<<",m_disabled="<<pln.m_disabled;
	out<<",m_bind_to_grid="<<pln.m_bind_to_grid;
	out<<",m_plane="<<pln.m_plane<<",span=("<<pln.m_uspan<<","<<pln.m_vspan<<")";
	out<<",num=("<<pln.m_unum<<","<<pln.m_vnum<<")";
	return out;
}
