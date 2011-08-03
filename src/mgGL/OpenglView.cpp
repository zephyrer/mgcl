/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// MGOpenGLView.cpp : インプリメンテーション ファイル
//
#include "MGCLStdAfx.h"
#include <utility>
#include "mg/DrawFunc.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Straight.h"
#include "mg/TrimmedCurve.h"
#include "mg/CSisect.h"
#include "mg/Group.h"
#include "mg/GelPositions.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/GLDrawFunc.h"
#include "mgGL/SysGLList.h"
#include "mgGL/glViewAttrib.h"
#include "mgGL/CommandDrawer.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"

using namespace mgGDL;
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace{ 
	size_t STRING_COUNT = 256; // is used to glGenLists and glDeleteLists.
}
#define INITIAL_SCALE 6.
#define INITIAL_BOX_SCALE INITIAL_SCALE*1.1
/////////////////////////////////////////////////////////////////////////////
// MGOpenGLView

/*double Bcolor2[3]={.6,.6,.6};//Background color.
double Gcolor2[3]={0.,0.,0.};//Object default color.
double Hcolor2[3]={1.,1.,0.};//Hilight color.*/
void MGOpenGLView::set_defalut_colors(){
	setBcolor(MGColor::get_instance(MGColor::Gray));
	setGcolor(MGColor::get_instance(MGColor::Black));
	setHcolor(MGColor::get_instance(MGColor::Yellow));
}

MGOpenGLView::MGOpenGLView(
	bool perspective	//indicates if the view is pespective or not.
):m_parent_glView(0), m_perspective(perspective),
m_object(0),m_line_density(1), 
m_center(0.,0.,0.),m_center_current(0.,0.,0.),m_scale(INITIAL_SCALE),
m_cx(0.),m_cy(0.),m_fovy(45.),m_diameter(1.),
m_hRC(0),m_hDC(0),m_smooth(.002), m_pick_aperture(5.),
m_eye_func(0),m_sitring_list_base(1),m_display_list(0){
//Smoothness is about 100 division of the screen's width curve:m_smooth(.01).
//Pick aperture is 5 pixels:m_pick_aperture(5.).

	set_defalut_colors();
	MGDrawFunc::set_functions(
		MGGLBeginLINE_STRIP,MGGLVertex,MGGLEnd,MGDrawPoint);
}

//Construct from MGglViewAttrib.
MGOpenGLView::MGOpenGLView(
	const MGglViewAttrib& glatr
):m_parent_glView(0), m_object(0),m_line_density(1),m_hRC(0),m_hDC(0),
m_eye_func(0),m_sitring_list_base(1),m_display_list(0),
m_perspective(glatr.m_perspective), m_center(glatr.m_center),
m_center_current(glatr.m_center),
m_scale(glatr.m_scale), m_cx(glatr.m_cx), m_cy(glatr.m_cy),
m_fovy(glatr.m_fovy),m_diameter(glatr.m_diameter),
m_cplane(glatr.m_cplane),m_eyeP(glatr.m_eyeP),m_up_vector(glatr.m_up_vector),
m_far(glatr.m_far),m_near(glatr.m_near),m_smooth(.01),
m_pick_aperture(5.){
//Smoothness is about 100 division of the screen's width curve:m_smooth(.01).
//Pick aperture is 5 pixels:m_pick_aperture(5.).

	set_defalut_colors();
	MGDrawFunc::set_functions(
		MGGLBeginLINE_STRIP,MGGLVertex,MGGLEnd,MGDrawPoint);

	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(glatr.m_currentMat);
}

MGOpenGLView::~MGOpenGLView(){
	if(!has_parent_OpenGLView()){
		make_RC_current();
		glDeleteLists(m_sitring_list_base, STRING_COUNT);
	}

	wglMakeCurrent(NULL,  NULL);
	wglDeleteContext(m_hRC);
	if( m_hDC )	::ReleaseDC(NULL,m_hDC);
	m_hDC = NULL;
}

//Copy the informations of glview2 into this.
//Data that is not copied from glview2 are:
//m_hRC(Rendering context), m_hDC(Device Context), m_object, m_display_list,
//m_eye_func, and m_sysgllist.
//m_hRC and m_hDC can be set using setDCRC(), m_object can be set by set_object().
//m_eye_func by set_eye_func each.
//m_display_list can be made using make_display_list().
//m_sysgllist must be made by invoking openGL's
//display list generation functions.
void MGOpenGLView::copy(const MGOpenGLView& glview2){
	m_cplane=glview2.m_cplane;

	// under construction
//	m_strlist=glview2.m_strlist;

	m_line_density=glview2.m_line_density;
	m_fovy=glview2.m_fovy;
	m_perspective=glview2.m_perspective;
	m_near=glview2.m_near;
	m_far=glview2.m_far;
	m_eyeP=glview2.m_eyeP;
	m_up_vector=glview2.m_up_vector;
	m_center_current=m_center=glview2.m_center;
	m_diameter=glview2.m_diameter;
	m_Bcolor=glview2.m_Bcolor;
	m_Gcolor=glview2.m_Gcolor;
	m_Hcolor=glview2.m_Hcolor;
	m_smooth=glview2.m_smooth;
	m_pick_aperture=glview2.m_pick_aperture;
	m_scale=glview2.m_scale;
	m_cx=glview2.m_cx;
	m_cy=glview2.m_cy;

	double currentMat[16];
	glview2.make_RC_current();
	glGetDoublev(GL_MODELVIEW_MATRIX,currentMat);
	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(currentMat);
}
void MGOpenGLView::copy(const MGglViewAttrib& glatr){
	m_center_current=m_center=glatr.m_center;
	m_cplane=glatr.m_cplane;//cout<<m_cplane<<endl;
	m_cx=glatr.m_cx;
	m_cy=glatr.m_cy;
	m_diameter=glatr.m_diameter;
	m_eyeP=glatr.m_eyeP;
	m_far=glatr.m_far;
	m_fovy=glatr.m_fovy;
	m_near=glatr.m_near;
	m_perspective=glatr.m_perspective;
	m_scale=glatr.m_scale;
	m_up_vector=glatr.m_up_vector;

	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(glatr.m_currentMat);
}

//Return display list name.
size_t MGOpenGLView::display_list()const{
	if(!has_parent_OpenGLView())
		return m_display_list;
	return get_parent_OpenGLView()->display_list();
}

//Draw the scene defined in this view including the current objects as hilighted.
void MGOpenGLView::drawScene(MGCL::VIEWMODE vmode, const MGPickObjects& pobjs){
	make_RC_current();
	if(has_display_list()){

	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT // | GL_POLYGON_BIT
				| GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT);//Push Attrib//////

	const float* Bcolr=m_Bcolor.color();
    glClearColor(Bcolr[0], Bcolr[1], Bcolr[2], Bcolr[3]);
	glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDepthFunc(GL_LEQUAL);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	set_frustum();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	set_model_matrix();//Set the model matrix.
	if(vmode==MGCL::SHADING || vmode==MGCL::WIRE_AND_SHADING)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);

//1. Construction plane drawing
	m_cplane.draw();
	
	MGOpenGLView* gltarget=this;
	if(m_parent_glView)
		gltarget=m_parent_glView;

//2. Target objects drawing
	m_Gcolor.exec();
	if(vmode==MGCL::SHADING || vmode==MGCL::WIRE_AND_SHADING)
		gltarget->drawShading();
	gltarget->drawWire(vmode,pobjs);
	
	glPopAttrib();//Pop Attrib////////

//3. the system generated display drawing.
	m_Gcolor.exec();
	gltarget->m_sysgllist.draw_list();

//4.draw command specific pictures.
	std::list<MGCommandDrawer*>::iterator i=gltarget->m_command_drawers.begin(),
		ie=gltarget->m_command_drawers.end();
	for(; i!=ie; i++){
		MGCommandDrawer* drawer=*i;
		if(drawer)
			drawer->draw();
	}

	glPopMatrix();
	::SwapBuffers(m_hDC);
	
	}//if(has_display_list())
	glFinish();
}

//Draw the scene defined in this view in shading mode.
void MGOpenGLView::drawShading(){
	mgGDL::setup_easy_shade();
	MGColor GcolorShading(m_Gcolor);
	GcolorShading+=float(.82);
	GcolorShading.exec();
	glCallList(m_display_list+2);
}

//Draw the scene defined in this view in wire mode, including the current objects as hilighted.
void MGOpenGLView::drawWire(MGCL::VIEWMODE vmode,const MGPickObjects& pobjs){
	mgGDL::setup_wire_drawing();

	glEnable(GL_POLYGON_OFFSET_FILL);// ポリゴンオフセットフィルを設定
	glPolygonOffset(1., 1.);

//1. the document drawing if wire drawing is necessary.
	m_Gcolor.exec();
	if(vmode==MGCL::WIRE || vmode==MGCL::WIRE_AND_SHADING)
		glCallList(m_display_list);

//2. current objects highlighting.
//**note** We are using (dlist_name()+1) as the highlighting display list name
//of the object gel().
	size_t nHL=pobjs.size();
	if(!nHL)
		return;

	double sl=span_length();
	int ld=line_density();
	glDisable(GL_DEPTH_TEST);
	for(size_t i=0; i<nHL; i++){
		pobjs[i].hilight_using_display_list(m_Hcolor,sl,ld);
	}
}

//get ModelView matrix of OpenGL.
void MGOpenGLView::get_model_matrix(
	double modelMat[16] //OpenGL's model matrix
)const{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	MGOpenGLView* this2=const_cast<MGOpenGLView*>(this);
	this2->set_model_matrix();//Set the model matrix.
	glGetDoublev(GL_MODELVIEW_MATRIX,modelMat);
	glPopMatrix();
}

//get projection matrix of OpenGL.
void MGOpenGLView::get_projection_matrix(
	double projMat[16]	//OpenGL's projection matrix
)const{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	set_frustum();
	glGetDoublev(GL_PROJECTION_MATRIX,projMat);
	glPopMatrix();
}

//get viewport of OpenGL.
void MGOpenGLView::get_viewport(
	int vp[4]	//OpenGL's projection matrix
)const{
	glGetIntegerv(GL_VIEWPORT,vp);
}

//Get the surface parameter value uv(u,v) where screen coordinate (sx,sy) is projected on.
//If no projection points are found, the nearest point to a perimeter of surf will
//be returned.
bool MGOpenGLView::get_surface_parameter_glv(
	const MGFSurface& surf,
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	MGPosition& uv	//surface parameter (u,v) where (sx,sy) is projected on will be returned.
)const{
	MGStraight sl;
	unproject_to_sl_glv(sx,sy,sl);
	MGCSisect_list csis=surf.isect(sl);
	if(csis.size()){
		uv=csis.front().param_surface();
		return true;
	}
	uv=surf.closest_on_boundary(sl);
	return false;
}

//Get the current viewing environment.
//This output can be used for the parameters of initialize_viewing_environment,
//and thus the current environment that were save can be restored.
void MGOpenGLView::get_viewing_environment(
	double& neard, double& fard,	//near and far data
	double& maxsize, MGPosition& center,	//Model's size(diameter) and the center
	double& scale, double& dx, double& dy,	//Current viewing's scale, and panning data(dx, dy).
	MGPosition& eye, MGVector& up_vector,//viewing eye position and view up vector
								//These parameters are used fro gluLookAt. see the explanation.
	double currentMat[16]	//OpenGL's ModelView current matrix.
)const{
	neard=m_near; fard=m_far;
	maxsize=m_diameter; center=m_center;
	scale=m_scale; dx=m_cx; dy=m_cy;
	eye=m_eyeP; up_vector=m_up_vector;
	make_RC_current();
	glGetDoublev(GL_MODELVIEW_MATRIX,currentMat);
}

//Test if this has a display list to draw.
bool MGOpenGLView::has_display_list()const{
	if(m_display_list) return true;
	if(!has_parent_OpenGLView()) return false;
	return get_parent_OpenGLView()->has_display_list();
}

//Transform to the home position.
void MGOpenGLView::home(){
	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	m_scale=INITIAL_SCALE;
	m_cx=0.; m_cy=0.;
	if(!m_perspective)
		return;

	MGMatrix mat;
	const MGPlane& cpl=m_cplane.plane();
	mat.set_xy_axis(cpl.u_deriv(), cpl.v_deriv());
	double glmat[16];
	mat.convert_to_glMatrix(glmat);
	glMultMatrixd(glmat);
}

//Initialize the viewing environment.
//Initialization is done by the parameter box or all the necessary environmental data.
//initialize_viewing_environment will invoke make_RC_current(), and so on return
//from initialize_viewing_environment, the rendering context is guaranteed to be current.
void MGOpenGLView::initialize_viewing_environment(
	const MGBox& box,
	int view_num
){
	m_center_current=m_center=box.mid();
	MGBox box2=box*INITIAL_BOX_SCALE;
	//std::cout<<box2.mid()<<box.mid()<<std::endl;///////
	m_diameter=box2.len();
	if(m_diameter<1.) m_diameter=1.;
	double diam10=m_diameter*10.;

	double tan_theta2;
	if(m_perspective)
		tan_theta2=tan(MGCL::degree_to_radian(m_fovy*.5))*2.;
	else
		tan_theta2=1.;
	m_near=diam10/tan_theta2;
	m_far=m_near+diam10;

	double eye=m_near+.5*diam10;
	if(m_eye_func)
		//Get the viewing environment data invoking m_eye_func.
		(*m_eye_func)(m_object,eye,m_eyeP,m_up_vector);
	else{
		m_eyeP=MGPosition(0.,0.,eye);
		m_up_vector=MGVector(0.,1.,0.);
	}
	m_cplane.set_grid_data(box2,view_num);
	home();
}

//Initialize the viewing environment to view objects projected onto a plane.
//This view is 2D view.
//center of the object is the origin of the plane.
void MGOpenGLView::initialize_viewing_environment(
	const MGPlane& plane,
	double height,	//height of the view.
					//If rotate=true, height is the height after rotated.
	bool rotate		//if true, 90 degree rotation will be performed
					//for the plane. (-v_deriv, u_deriv) will be (x,y).
					//If false, (u_deriv, v_deriv) will be (x,y).
){
	m_perspective=false;
	m_center_current=m_center=plane.root_point();
	m_diameter=height*INITIAL_BOX_SCALE;
	if(m_diameter<1.) m_diameter=1.;
	double diam10=m_diameter*10.;

	m_near=diam10;
	m_far=diam10+diam10;
	home();
	m_eyeP=MGPosition(0.,0.,m_near+diam10*.5);
	m_up_vector=MGPosition(0.,1.,0.);
	MGMatrix mat;
	if(rotate){
		mat.set_xy_axis(-plane.v_deriv(),plane.u_deriv());
	}else{
		mat.set_xy_axis(plane.u_deriv(), plane.v_deriv());
	}
	double glmat[16];
	mat.convert_to_glMatrix(glmat);
	glMultMatrixd(glmat);
}

void MGOpenGLView::initialize_viewing_environment(
	double neard, double fard,	//near and far data
	double maxsize, const MGPosition& center,	//Model's size(diameter) and the center
	double scale, double dx, double dy,	//Current viewing's scale, and panning data(dx, dy).
	const MGPosition& eye, const MGVector& up_vector,//viewing eye position and view up vector
								//These parameters are used fro gluLookAt. see the explanation.
	const double currentMat[16]	//OpenGL's ModelView current matrix.
){
	m_near=neard; m_far=fard;
	m_diameter=maxsize; m_center_current=m_center=center;
	m_scale=scale; m_cx=dx; m_cy=dy;
	m_eyeP=eye; m_up_vector=up_vector;
	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(currentMat);
}
void MGOpenGLView::initialize_viewing_environment(
	const MGglViewAttrib& viewA
){
	copy(viewA);
}

//locate the screen coordinate (x,y) in the  3D world coordinate.
//(x, y)'s origin is (left, bottom) of the screen.
MGPosition MGOpenGLView::locate_glv(int x, int y, MGPosition* uv)const{
	MGStraight sl;
	unproject_to_sl_glv(x,y,sl);
	const MGConstructionPlane& pl=construction_plane();
	MGPosition xyz,uv2;
	if(pl.enabled()) xyz=pl.locate(sl,uv2);
	else xyz=sl.root_point();

	if(uv) *uv=uv2;
	return xyz;
}

//construct the construction plane along with its display list.
void MGOpenGLView::make_construction_plane(
	const MGBox& bx,
	int view_num	//Standard view number:
	//1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	//0: non standard view.
){
	m_cplane.set_grid_data(bx,view_num);
}

//construct the construction plane along with its display list.
void MGOpenGLView::make_construction_plane(
	const MGPosition& mid,	//center of the construction plane.
	const MGVector& uderi,	//u-axis vector of the construction plane.
	const MGVector& vderi,	//v-axis vector of the construction plane.
	double uspan,			//span length between the lines along u-axis.
	double vspan,			//span length between the lines along v-axis.
	size_t ulnum,			//number of lines to draw along u-axis.
	size_t vlnum			//number of lines to draw along v-axis.
){
	m_cplane.set_grid_data(MGPlane(uderi,vderi,mid), uspan,vspan,ulnum,vlnum);
}

//Make openGL display list in this glview.
size_t MGOpenGLView::make_display_list(const MGGroup& grp){
//Generate OpenGL display list.
	if(has_parent_OpenGLView())
		return get_parent_OpenGLView()->m_display_list;
	UseFontBitmaps();
	m_display_list=grp.make_display_list(span_length(),line_density());
	return m_display_list;
}

//Function's return value is the number of hit objects.
int MGOpenGLView::pick_to_select_buf(
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	double aperturex, double aperturey,
	size_t display_list,	//display list that includes pick objects.
	size_t buf_size, GLuint* selectBuf,
	GLint* vport,		//view port will be output for the pick operation.
	GLdouble* modelMatrix,//model matrix will be output for the pick operation.
	GLdouble* projMatrix	//projection matrix will be output for the pick operation.
){
	glSelectBuffer(buf_size,selectBuf);
	GLint viewport[4];
	GLint* vp=vport; if(!vp) vp=viewport;
	glGetIntegerv(GL_VIEWPORT,vp);
	glRenderMode(GL_SELECT);

	glInitNames();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPickMatrix((GLdouble)sx,(GLdouble)sy,aperturex,aperturey,vp);
	set_frustum();
	if(projMatrix)
		glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
		set_model_matrix();//set_model_matrix
		//Draw objects to pick.
		glCallList(display_list);
		if(modelMatrix)
			glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	int objnum=	glRenderMode(GL_RENDER);
	if(objnum<0){
		//This occurs when num of picked objected exceeds BUFSIZE.
		//In this case, objnum=-1 was returned.
		size_t bufCount=0;
		objnum=-1;
		while(bufCount<=buf_size){
			objnum++;
			bufCount+=selectBuf[bufCount];
			bufCount+=3;
		};
	}
	return objnum;
}

#define BUFSIZE 1024
//Pick objects in the display list generated by make_display_list.
//Function's return value is MGPickObject vector in m_CurrentObjects member data.
//All the objects which were inside the pick aperture will be output.
//This data can be accessed using current_object(), or current_PickObject().
//pick will invoke make_RC_current();
MGPickObjects MGOpenGLView::pick_glv(
	MGCL::VIEWMODE vmode,//View mode of the document.
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	double aperturex, double aperturey,
	const MGAbstractGels& objtypes,
	bool get_param
){
	MGPickObjects pobjs;
	if(!has_display_list())
		return pobjs;

	GLuint selectBuf[BUFSIZE];
	if(aperturex<=0.)
		aperturex=pick_aperture();
	if(aperturey<=0.)
		aperturey=pick_aperture();

	make_RC_current();
	size_t dlname=display_list();
	if(vmode==MGCL::SHADING || vmode==MGCL::WIRE_AND_SHADING)
		dlname+=2;
	int objnum=pick_to_select_buf(
				sx,sy,aperturex,aperturey,dlname,BUFSIZE,selectBuf);

	MGStraight sl;
	unproject_to_sl_glv(sx,sy,sl);
	size_t idBuf=0;
	int npick;
	for(int i=0; i<objnum; i++, idBuf+=npick){
		npick=selectBuf[idBuf];
		if(!npick)
			continue;

		idBuf+=3;
		int npickm1=npick-1;
		MGPickObject pobj;
		for(int j=0; j<npickm1; j++){
			MGGel* grp=(MGGel*)(selectBuf[idBuf+j]);
			pobj.set_group(static_cast<MGGroup*>(grp));
		}

		MGGel* gl=(MGGel*)(selectBuf[idBuf+npickm1]);
		MGCellNB* cell=dynamic_cast<MGCellNB*>(gl);
		MGComplex* cmplx=0;
		if(cell)
			cmplx=cell->parent_complex();
		
		if(cmplx){
			if(cmplx->type_is(objtypes)){
				if(!pobjs.includes(cmplx)){
					pobj.set_gel(cmplx);
					pobjs.push_back(pobj);
				}
			}
			continue;
		}else if(!gl->type_is(objtypes))
			continue;

		pobj.set_gel(gl);
		if(get_param){
			const MGObject* obj=gl->includes_object();
			if(obj)
				pobj.set_parameter(obj->pick_closest(sl));
		}
		pobjs.push_back(pobj);
		
	}

	return pobjs;
}

//Pick a perimeter of the surface surf. That is, obtain the perimeter number
//that passes input (sx,sy) when drawn in the current view matrix.
//Function's return value is perimeter number picked.
//When no perimeters are picked, -1 will be returned.
int MGOpenGLView::pick_perimeter_glv(
	const MGSurface& surf,
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv,	//surface parameter (u,v) nearest to (sx,sy) will be returned.
	double aperturex,//specifies pick aperture of x and y.
	double aperturey//When <=0. value is specified, default value(the value
			//obtained by pick_aperture() will be used.
){
	int perimeter=-1;
	if(aperturex<=0.) aperturex=pick_aperture();
	if(aperturey<=0.) aperturey=pick_aperture();

	GLuint selectBuf[BUFSIZE];
	glSelectBuffer(BUFSIZE,selectBuf);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);

	glInitNames();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	gluPickMatrix((GLdouble)sx,(GLdouble)sy,aperturex,aperturey,viewport);
	set_frustum();
	GLdouble projMatrix[16];
	glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

	size_t nperi=surf.perimeter_num();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
		set_model_matrix();//set_model_matrix
		//Draw objects to pick.
		for(size_t i=0; i<nperi; i++){
			MGCurve* peri=surf.perimeter_curve(i);
			glRenderMode(GL_SELECT);
			glPushName(i);
			peri->drawWire(span_length(),line_density());
			
			glPopName();
			int objnum=	glRenderMode(GL_RENDER);
			if(objnum>0){
				perimeter=int(i);
				
				if(uv){ // if parameter is needed
					GLdouble modelMatrix[16];
					glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
					double t;
					get_near_position(peri,sx,sy,modelMatrix,projMatrix,viewport,t);
					*uv=surf.perimeter_uv(perimeter,t);
				}

				delete peri;
				break;
			}
			delete peri;
		}
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	return perimeter;
}

//Pick an edge of the face f. That is, obtain the edge number
//that passes input (sx,sy) when drawn in the current view matrix.
//Function's return value is the edge pointer picked.
//When no edges are picked, null will be returned.
const MGEdge* MGOpenGLView::pick_edge_glv(
	const MGFace& f,
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv,	//surface parameter (u,v) nearest to (sx,sy) will be returned.
	double aperturex,//specifies pick aperture of x and y.
	double aperturey//When <=0. value is specified, default value(the value
			//obtained by pick_aperture() will be used.
){
	const MGEdge* edge=0;
	if(aperturex<=0.) aperturex=pick_aperture();
	if(aperturey<=0.) aperturey=pick_aperture();

	GLuint selectBuf[BUFSIZE];
	glSelectBuffer(BUFSIZE,selectBuf);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);

	glInitNames();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPickMatrix((GLdouble)sx,(GLdouble)sy,aperturex,aperturey,viewport);
	set_frustum();
	GLdouble projMatrix[16];
	glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
		set_model_matrix();//set_model_matrix
		//Draw objects to pick.
		size_t nloop=f.number_of_loops();
		for(size_t i=0; i<nloop; i++){
			const MGLoop& li=*(f.loop(i));
			size_t nedge=li.number_of_edges();
			for(size_t j=0; j<nedge; j++){
				const MGEdge* edge2=li.edge(j);
				MGEdge& be=*(edge2->make_binder_with_curve());
				MGTrimmedCurve cij=be.trimmed_curve();
				glRenderMode(GL_SELECT);
				glPushName(j);
				cij.drawWire(span_length(),line_density());
				glPopName();
				int objnum=	glRenderMode(GL_RENDER);
				if(objnum>0){
					edge=edge2;
					GLdouble modelMatrix[16];
					glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
					double t;
					get_near_position(&cij,sx,sy,modelMatrix,projMatrix,viewport,t);
					if(uv) *uv=edge2->eval(be.param_pcell(t));
					break;
				}
			}
			if(edge) break;
		}
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	return edge;
}

//Determine if screen coordinate (sx,sy) is closer to the start point or to the end
//of the curve curve.
//Functin's return value is 0: if start point, 1: if end point.
int MGOpenGLView::pick_start_end_glv(
	const MGCurve& curve,
	int sx, int sy	//Screen coordinates. (left, bottom) is (0,0).
){
	MGStraight sl;
	unproject_to_sl_glv(sx,sy,sl);
	MGPosition P0=curve.start_point(), P1=curve.end_point();
	if(sl.distance(P0)<=sl.distance(P1)) return 0;
	return 1;
}

//Project world coordinates to OpenGL's screen coordinates.
//If modelMat, projMat, or vp is not input, project will ask OpenGL to get them.
//Generally, users of project are recommended to get modelMat, projlMat, or
//vp, and input them to project if continuous multiple use of project will take place.
//If one of modelMat, projlMat, or vp is not input, make_RC_current() must be invoked
//before use of project.
void MGOpenGLView::project(
	const MGPosition& world, MGPosition& screen,
	const double* modelMat,	//OpenGL's model matrix
	const double* projlMat,	//OpenGL's projection matrix
	const int* vp			//OpenGL's viewport
) const{
	const double* model=modelMat;
	const double* proj=projlMat;
	const int* vp2=vp;
	GLdouble modelM[16], projM[16];
	GLint vp3[4];
	if(!proj){
		proj=projM;
		get_projection_matrix(projM);
	}
	if(!model){
		model=modelM;
		get_model_matrix(modelM);
	}
	if(!vp2){
		vp2=vp3;
		glGetIntegerv(GL_VIEWPORT,vp3);
	}

	screen.resize(3);
	double* x=screen.data();
	int ret=gluProject(world[0],world[1],world[2],model,proj,vp2,x,x+1,x+2);
}

//Rotate the current view by the angle along the vector(x,y,z),
//performs a counterclockwise rotation of angle angle about
//the vector from the origin through the point (x, y, z).
void MGOpenGLView::rotate(double angle, double x, double y, double z){
	make_RC_current();
	glMatrixMode(GL_MODELVIEW);
	double saveMat[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,saveMat);
	glLoadIdentity();

	// 拡大・縮小
	glTranslated(m_cx,m_cy,0.);
	glRotated(angle, x,y,z);//Rotate around scene origin.
	glTranslated(-m_cx,-m_cy,0.);
	glMultMatrixd(saveMat);
}

void MGOpenGLView::set_center_current(int x, int y){
	GLdouble modelMatrix[16], projMatrix[16];
	get_projection_matrix(projMatrix);
	get_model_matrix(modelMatrix);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	MGStraight sl;
	unproject_to_sl(x,y,modelMatrix,projMatrix,viewport,sl);
	m_center_current=sl.eval(sl.param(m_center));
}

//Set background color;
void MGOpenGLView::setBcolor(const MGColor& color){
	m_Bcolor=color;
}

//Set default object color;
void MGOpenGLView::setGcolor(const MGColor& color){
	m_Gcolor=color;
}

//Set hilight color;
void MGOpenGLView::setHcolor(const MGColor& color){
	m_Hcolor=color;
}

//Get redering context.
void MGOpenGLView::setDCRC(HDC dc, HGLRC rc){
	if(m_hRC) wglMakeCurrent(m_hDC, NULL);
	m_hDC=dc;
	m_hRC=rc;
	make_RC_current();
}

//Set Frustum of OpenGL.
void MGOpenGLView::set_frustum()const{
	double height2=view_volume_height()*.5;
	GLint vp[4]; glGetIntegerv(GL_VIEWPORT,vp);
	double wide2=height2*double(vp[2])/double(vp[3]);
	double left=m_cx-wide2, right=m_cx+wide2;
	double bottom=m_cy-height2, top=m_cy+height2;

//	double znear=view_volume_near()/2., zfar=view_volume_far()*2.;
	double znear=view_volume_near(), zfar=view_volume_far();
	if(m_perspective)
		glFrustum(left,right,bottom,top,znear,zfar);
	else
		glOrtho(left,right,bottom,top,znear,zfar);
}

//Set Model_View and Projection matrix of OpenGL.
void MGOpenGLView::set_model_matrix(){
	double currentMat[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,currentMat);
	glLoadIdentity();

	const MGPosition& eye=eye_position();
	const MGVector& up=view_up_vector();
	gluLookAt(eye[0],eye[1],eye[2],0.,0.,0., up[0],up[1],up[2]);
	m_sysgllist.pre_transform_process();
	glMultMatrixd(currentMat);
	const MGPosition& cntr=center();
	glTranslated(-cntr[0], -cntr[1], -cntr[2]);
}

//Set the parent MGOpenGLView.
void MGOpenGLView::set_parent_OpenGLView(MGOpenGLView* parent){
	if(parent){
		ASSERT(m_parent_glView==0);
		wglShareLists(parent->getRC(),getRC());
	}
	m_parent_glView=parent;
}

//Convert the screen coordinate (sx, sy) to world coordinate (wx, wy) on the 
//view plane.
void MGOpenGLView::screen_to_world(
	int wh[2],	//width(wh[0]) and height(wh[1]) of the screen.
	double sx,double sy, double& wx, double& wy
)const{
	double sheight=double(wh[1]), swidth=double(wh[0]);
	double wsratio=view_volume_height()/sheight;
	wx=m_cx+wsratio*(sx-swidth*.5);
	wy=m_cy+wsratio*(sy-sheight*.5);
}

//Translate the current view by (dx, dy).
void MGOpenGLView::translate(double dx, double dy){
	m_cx-=dx/m_scale; m_cy-=dy/m_scale;
}

//Translate the current view by (dx, dy) without current scale.
void MGOpenGLView::translate_without_scale(double dx, double dy){
	m_cx-=dx; m_cy-=dy;
}

//Translate and scale the current view.
//(x0, y0) to (x1,y1) is the rectangle of screen coordinate whose origin is
//(left,bottom).
void MGOpenGLView::pan_zoom(int x0, int y0, int x1, int y1){
	make_RC_current();
	int dx=x1-x0, dy=y1-y0;
	if(dx==0 && dy==0) return;

	if(dx<0) dx*=-1;if(dy<0) dy*=-1;
	double sdy=double(dy);

	GLint vp[4]; glGetIntegerv(GL_VIEWPORT,vp);
	double sheight=double(vp[3]), swidth=double(vp[2]);
	if(dx){
		double sdx=double(dx);
		double aspect=sheight/swidth;
		if(sdy/sdx < aspect) sdy=sdx*aspect;
	}

	double sxm=(double(x0)+double(x1))*.5;
	double sym=(double(y0)+double(y1))*.5;
	screen_to_world(vp+2,sxm,sym,m_cx,m_cy);
	double wsratio=view_volume_height()/sheight;
	m_scale=diameter()/(wsratio*sdy);
}

//Translate and scale the current view.
//box is the world coordinate's box cube.
void MGOpenGLView::pan_zoom(const MGBox& box){
	make_RC_current();
	MGPosition wc=box.mid();
	MGPosition sc;
	double mmat[16]; get_model_matrix(mmat);
	int vp[4]; get_viewport(vp);
	project(wc,sc,mmat,0,vp);
	screen_to_world(vp+2,sc[0],sc[1],m_cx,m_cy);
	double p=mmat[2]*wc[0]+mmat[6]*wc[1]+mmat[10]*wc[2]+mmat[14];
	p*=-1.;
	double len=view_volume_near()*box.len()/p;
	m_scale=diameter()/len;
}

//Convert the windows screen coordinate (x,y) to MGCL's straight line.
//and get the intersection of the straight line and the construction plane.
//The origin of the screen coordinate is left, bottom. Not left, top.
void MGOpenGLView::unproject(
	int x, int y,	//screen coordinate whose origin is (left, bottom).
	MGStraight& sl,	//The straight line of (x,y) will be returnred.
	MGCSisect& is	//the intersectio of the sl and the construction plane
					//will be returned.
)const{
	unproject_to_sl_glv(x,y,sl);
	if(m_cplane.valid()) sl.relation(m_cplane.plane(), is);
	else is=MGCSisect(sl.root_point(), 0.,MGPosition(0.,0.));
}

//Convert the windows screen coordinate (x,y) to MGCL's straight line.
//The origin of the screen coordinate is left, bottom. Not left, top.
void MGOpenGLView::unproject_to_sl_glv(int x, int y, MGStraight& sl)const{
	make_RC_current();
	GLdouble modelMatrix[16], projMatrix[16];
	get_projection_matrix(projMatrix);
	get_model_matrix(modelMatrix);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	::unproject_to_sl(x,y,modelMatrix,projMatrix,viewport,sl);
}

//Convert the windows screen coordinate (x,y) to MGCL's straight line.
//The origin of the screen coordinate is left, bottom. Not left, top.
void unproject_to_sl(
	int x, int y,
	const double modelMatrix[16],	//OpneGL's model matrix.
	const double projMatrix[16],	//OpneGL's projection matrix.
	const int viewport[4],			//OpneGL's viewport.
	MGStraight& sl					//straight will be output.
){
	MGPosition origin(3), point(3);
	double *uvw0=origin.data(),*uvw1=point.data();
	double xd=double(x), yd=double(y);
	int ret=gluUnProject(xd, yd, 0.,
			modelMatrix,projMatrix,viewport,uvw0,uvw0+1,uvw0+2);
	ret=gluUnProject(xd, yd, 100.,
			modelMatrix,projMatrix,viewport,uvw1,uvw1+1,uvw1+2);
	sl=MGStraight(MGSTRAIGHT_UNLIMIT, point-origin,origin);
}

void get_near_position(
	const MGCurve* crv, 
	int sx,int sy,	//screen coordinates whose origin is (left, bottom).
	const double modelMatrix[16],	//OpneGL's model matrix.
	const double projMatrix[16],	//OpneGL's projection matrix.
	const int viewport[4],			//OpneGL's viewport.
	double& t	//parameter value of the curve crv near to (sx,sy) will be returned.
){
	MGStraight sl;
	::unproject_to_sl(sx, sy, modelMatrix,projMatrix,viewport,sl);
	MGCurve* crv2=crv->clone();//std::cout<<"1:"<<(*crv2)<<std::endl;
	*crv2-=sl.root_point();
	MGMatrix mat; mat.set_axis(sl.direction(),2);
	*crv2*=mat;// std::cout<<"2:"<<(*crv2)<<std::endl;

	double dist;
	t=crv2->closest2D(MGDefault::origin_2D(),dist);

	delete crv2;
}

// grid snap enabled
void MGOpenGLView::enable_grid_snap(bool bEnabled){
	MGConstructionPlane& cplane=construction_plane();
	if(!cplane.valid()) return;

	if(cplane.disabled()) cplane.set_enable();
	if(bEnabled)
		cplane.set_bind_to_grid_enable();
	else
		cplane.set_bind_to_grid_disable();
}

// grid snap is enabled
bool MGOpenGLView::is_enabled_grid_snap()const{
	if(construction_plane().disabled()) return false;
	return construction_plane().is_bind_to_grid();
}

//Display string at the position pos with the color
void MGOpenGLView::draw_string(
	const CString& str,	//String to display
	const MGPosition& pos,//position for the string to display at.
	const MGColor* color//color. When color=null, white color will be employed.
){
	if(str.IsEmpty()) return;

	const MGColor* dcolor=&(MGColor::get_instance(MGColor::White));
	if(color)
		dcolor=color;
	dcolor->exec();
	glRasterPos3dv(pos.data());
	glListBase(m_sitring_list_base);
	glCallLists(str.GetLength(), GL_UNSIGNED_BYTE, (LPCTSTR)str);
}

//Invoke wglUseFontBitmaps(getHDC(),0,STRING_COUNT,m_sitring_list_base);
bool MGOpenGLView::UseFontBitmaps(){
	// make the system font the device context's selected font 
	HDC hdc=getHDC();
	//SelectObject (hdc, GetStockObject (SYSTEM_FONT));  
	BOOL result=wglUseFontBitmaps(hdc,0,STRING_COUNT,m_sitring_list_base);
	//DWORD error=GetLastError();
	return result==TRUE;
}
