/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Context.h"
#include "mgGL/OpenGLView.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////
// MGContext defines the attributes of a document.

MGContext::MGContext():m_line_density(1),m_smooth(float(.01)),
m_pick_aperture(float(5.)),m_appearance(0),m_gel(0),
//m_Bcolor(MGColor::get_instance(MGColor::DarkGray)),//Background color.
m_Bcolor(.5,.5,.5),//Background color.
m_Gcolor(MGColor::get_instance(MGColor::Black)),//Object lines color.
m_Hcolor(MGColor::get_instance(MGColor::Yellow))//Object highlight color.
{
	for(int i=0; i<4; i++) m_views[i]=0;

	m_tolerance[0]=0.5E-3;			//wc_zero
	m_tolerance[1]=1.0E-6;			//rc_zero
	m_tolerance[2]=1.0E-20;			//mach_zero
	m_tolerance[3]=m_tolerance[0];	//line_zero;
	m_tolerance[4]=.0025;			//angle_zero;
	m_tolerance[5]=5.0E+2;			//max_knot_ratio;
}

MGContext::MGContext(
	const MGSnapAttrib& snap_attrib,//Snap data
	
	int line_density,		//line density for a surface to draw in wire mode.
	const MGColor& Bcolor,	//Background color.
	const MGColor& Gcolor,	//Object lines color.
	const MGColor& Hcolor,	//Object highlight color.
	float smooth,	//Smoothness of the curves to draw.
		// 1/smooth is the division number of a curve whose length is the window width.
		// When smooth becomes small, smoothness increases.	
	float pick_aperture,//Pick aperture. Number of pixels to allow picking.
	MGglViewAttrib* view[4],//Newed objects of MGglViewAttrib, may be null.
		//The ownership will be transfered to this MGContext.
		//[0]=for 3D perspective view
		//[1]=for (x,y) 2D view
		//[2]=for (y,z) 2D view
		//[3]=for (z,x) 2D view

	const double torelance[6],//MGCL Tolerance.
		//[0]=wc_zero;
		//[1]=rc_zero;
		//[2]=mach_zero;
		//[3]=line_zero;
		//[4]=angle_zero;
		//[5]=max_knot_ratio;
	const mgTLInputParam& tessellate_param,//tessellation parameter.
	MGAppearance* appearance//must be a newed object of MGAppearance.
		//the ownership will be transfered to this MGContext.
):m_snap_attrib(snap_attrib),m_line_density(line_density),m_smooth(smooth),
m_pick_aperture(pick_aperture),m_tessellate_param(tessellate_param),
m_appearance(appearance),m_gel(0),
m_Bcolor(Bcolor),	//Background color.
m_Gcolor(Gcolor),	//Object lines color.
m_Hcolor(Hcolor)	//Object highlight color.
{
	size_t i;
	for(i=0; i<4; i++) m_views[i]=view[i];
	for(i=0; i<6; i++) m_tolerance[i]=torelance[i];
}

//copy constructor
MGContext::MGContext(
	const MGContext& context
):m_snap_attrib(context.m_snap_attrib),m_line_density(context.m_line_density),
m_smooth(context.m_smooth),m_pick_aperture(context.m_pick_aperture),
m_tessellate_param(context.m_tessellate_param)
,m_appearance(0),m_gel(0),
m_Bcolor(context.m_Bcolor),	//Background color.
m_Gcolor(context.m_Gcolor),	//Object lines color.
m_Hcolor(context.m_Hcolor)	//Object highlight color.
{
	if(context.m_appearance)
		m_appearance=context.m_appearance->clone();
	if(context.m_gel)
		m_gel=context.m_gel->clone();

	size_t i;
	for(i=0; i<4; i++){
		m_views[i]=0;
		if(context.m_views[i])
			m_views[i]=new MGglViewAttrib(*(context.m_views[i]));
	}
	for(i=0; i<6; i++)
		m_tolerance[i]=context.m_tolerance[i];
}

//Assignment operator
MGContext& MGContext::operator=(const MGContext& context){
	if(this==&context)
		return *this;

	m_snap_attrib=context.m_snap_attrib;
	m_line_density=context.m_line_density;
	m_smooth=context.m_smooth;
	m_pick_aperture=context.m_pick_aperture;
	m_tessellate_param=context.m_tessellate_param;

	delete m_appearance; m_appearance=0;
	if(context.m_appearance)
		m_appearance=context.m_appearance->clone();
	if(context.m_gel)
		m_gel=context.m_gel->clone();

	delete m_gel; m_gel=0;
	if(context.m_gel)
		m_gel=context.m_gel->clone();

	m_Bcolor=context.m_Bcolor;	//Background color.
	m_Gcolor=context.m_Gcolor;	//Object lines color.
	m_Hcolor=context.m_Hcolor;	//Object highlight color.

	size_t i;
	for(i=0; i<4; i++){
		delete m_views[i]; m_views[i]=0;
		if(context.m_views[i]) m_views[i]=new MGglViewAttrib(*(context.m_views[i]));
	}
	for(i=0; i<6; i++) m_tolerance[i]=context.m_tolerance[i];
	return *this;
}
MGContext& MGContext::operator=(const MGGel& gel2){
	const MGContext* gel2_is_this=dynamic_cast<const MGContext*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGContext::operator<(const MGContext& gel2)const{
	return this<&gel2;
}
bool MGContext::operator<(const MGGel& gel2)const{
	const MGContext* gel2_is_this=dynamic_cast<const MGContext*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//////////// Destructor.////////

MGContext::~MGContext(){
	delete m_appearance;
	for(size_t i=0; i<4; i++) delete m_views[i];
	delete m_gel;
}

void MGContext::set_Bcolor(const MGColor& Bcolor){
	m_Bcolor=Bcolor;	//Background color.
}
void MGContext::set_Gcolor(const MGColor& Gcolor){
	m_Gcolor=Gcolor;	//Background color.
}
void MGContext::set_Hcolor(const MGColor& Hcolor){
	m_Hcolor=Hcolor;	//Background color.
}

//Set the viewing context of MGOpenGLView view to the document.
void MGContext::set_view_context(
	int view_num,	//Standard view number of the view.
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
	const MGOpenGLView& view
){
	set_line_density(view.line_density());
	set_smooth(view.smooth());
	set_pick_aperture(view.pick_aperture());
	set_Bcolor(view.Bcolor());
	set_Gcolor(view.Gcolor());
	set_Hcolor(view.Hcolor());
	if(view_num<=0 || view_num>4) return;
	set_view(new MGglViewAttrib(view), view_num);
}

//Set the view data of the view view_num
void MGContext::set_view(
	MGglViewAttrib* view,//must be a newed object of MGglViewAttrib.
		//the ownership will be transfered to this MGContext.
	int view_num	//Standard view number:
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
){
	if(view_num<=0 || view_num>4) return;
	view_num--;
	delete m_views[view_num];
	m_views[view_num]=view;
};

////////get MGglViewAttrib attrib data.////////
//When view_num is not the number of a standard view, null will be returned.
const MGglViewAttrib*  MGContext::view(
	int view_num	//Standard view number of the view.
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
)const{
	if(view_num<=0 || view_num>4) return 0;
	view_num--;
	return m_views[view_num];
}
MGglViewAttrib*  MGContext::view(
	int view_num	//Standard view number of the view.
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
){
	if(view_num<=0 || view_num>4) return 0;
	view_num--;
	return m_views[view_num];
}

//Execution of MGCL tolerance.
void MGContext::exec_tolerance()const{
	MGTolerance::set_wc_zero(m_tolerance[0]);//
	MGTolerance::set_rc_zero(m_tolerance[1]);//
	MGTolerance::set_mach_zero(m_tolerance[2]);//
	MGTolerance::set_line_zero(m_tolerance[3]);//
	MGTolerance::set_angle_zero(m_tolerance[4]);//
	MGTolerance::set_max_knot_ratio(m_tolerance[5]);//
}

//Tolerance
void MGContext::set_tolerance(
	double wc_zero,
	double rc_zero,
	double mach_zero,
	double line_zero,
	double angle_zero,
	double max_knot_ratio
){
	m_tolerance[0]=wc_zero;
	m_tolerance[1]=rc_zero;
	m_tolerance[2]=mach_zero;
	m_tolerance[3]=line_zero;
	m_tolerance[4]=angle_zero;
	m_tolerance[5]=max_knot_ratio;
}

//Appearance data.
//set appearance.
void MGContext::set_appearance(
	MGAppearance* appearance	//appearance must be a newed object. The ownership will
								//be transfered to this MGContext.
){
	delete m_appearance;
	m_appearance=appearance;
}

//Execution of drawing OpenGL attributes.
//Valid OpenGL rendering context must be made current.
void MGContext::exec_draw_attributes()const{
	glDisable(GL_LIGHTING);
	if(!m_appearance) return;
	m_appearance->drawAttrib();
}

//Execution of rendering OpenGL attributes.
//Valid OpenGL rendering context must be made current.
void MGContext::exec_render_attributes()const{
	if(!m_appearance) return;
	m_appearance->render();
}

// Output virtual function.
std::ostream& MGContext::out(std::ostream& ostrm) const{
	const float* Bcolor=m_Bcolor.color();
	const float* Gcolor=m_Gcolor.color();
	const float* Hcolor=m_Hcolor.color();
	ostrm<<std::endl<<"MGContext="<<this<<","<<m_snap_attrib<<std::endl;
	ostrm<<",Bcolor["<<Bcolor[0]<<","<<Bcolor[1]<<","<<Bcolor[2]<<"]";
	ostrm<<",Gcolor["<<Gcolor[0]<<","<<Gcolor[1]<<","<<Gcolor[2]<<"]";
	ostrm<<",Hcolor["<<Hcolor[0]<<","<<Hcolor[1]<<","<<Hcolor[2]<<"]";
	ostrm<<",smooth="<<m_smooth;
	ostrm<<",pick_aperture="<<m_pick_aperture;

	size_t i;
	for(i=0; i<4; i++){
		ostrm<<",views["<<i<<"]=";
		if(m_views[i]){
			ostrm<<(*(m_views[i]));
		}else{
			ostrm<<"NULL";
		}
	}
	
	ostrm<<",wc_zero="<<m_tolerance[0];
	ostrm<<",rc_zero="<<m_tolerance[1];
	ostrm<<",mach_zero="<<m_tolerance[2];
	ostrm<<",line_zero="<<m_tolerance[3];
	ostrm<<",angle_zero="<<m_tolerance[4];
	ostrm<<",max_knot_ratio="<<m_tolerance[5];

	ostrm<<",tessellate_param="<<m_tessellate_param;

	ostrm<<",appearance=";
	if(m_appearance){
		ostrm<<(*m_appearance);
	}else{
		ostrm<<"NULL";
	}

	ostrm<<",m_gel=";
	if(m_gel){
		ostrm<<(*m_gel);
	}else{
		ostrm<<"NULL";
	}

	return ostrm;
}

//Write all member data.
void MGContext::WriteMembers(MGOfstream& buf)const{
	size_t i;

	buf<<m_snap_attrib;
	buf<<m_line_density;
	const float* Bcolor=m_Bcolor.color();
	for(i=0; i<3; i++) buf<<Bcolor[i];
	const float* Gcolor=m_Gcolor.color();
	for(i=0; i<3; i++) buf<<Gcolor[i];
	const float* Hcolor=m_Hcolor.color();
	for(i=0; i<3; i++) buf<<Hcolor[i];
	buf<<m_smooth;
	buf<<m_pick_aperture;
	for(i=0; i<4; i++){
		if(m_views[i]){
			buf<<0xffffffffL;
			//The header that indicates an MGglViewAttrib object is followed.
			buf<<(*(m_views[i]));
		}else{
			//When null pointer.
			buf<<0x00000000L;
		}
	}
	for(i=0; i<6; i++) buf<<m_tolerance[i];
	buf<<m_tessellate_param;

	buf.WritePointer(m_appearance);
	buf.WritePointer(m_gel);
}

//Write all member data
void MGContext::ReadMembers(MGIfstream& buf){
	size_t i;

	buf>>m_snap_attrib;
	buf>>m_line_density;
	float* Bcolor=m_Bcolor.color();
	for(i=0; i<3; i++) buf>>Bcolor[i]; Bcolor[3]=1.;
	float* Gcolor=m_Gcolor.color();
	for(i=0; i<3; i++) buf>>Gcolor[i]; Gcolor[3]=1.;
	float* Hcolor=m_Hcolor.color();
	for(i=0; i<3; i++) buf>>Hcolor[i]; Hcolor[3]=1.;
	buf>>m_smooth;
	buf>>m_pick_aperture;
	int addr;
	for(i=0; i<4; i++){
		buf>>addr;
		if(addr){
			//The header that indicates an MGglViewAttrib object is followed.
			m_views[i]=new MGglViewAttrib;
			buf>>(*(m_views[i]));
		}else{
			//When null pointer.
			m_views[i]=0;
		}
	}
	for(i=0; i<6; i++) buf>>m_tolerance[i];
	buf>>m_tessellate_param;

	m_appearance=static_cast<MGAppearance*>(buf.ReadPointer());
	m_gel=buf.ReadPointer();

}
