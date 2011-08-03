/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgOpenGLView_HH_
#define _mgOpenGLView_HH_

class MGBox;
class MGGroup;
class MGObject;
class MGGLAttrib;
class MGAttribedGel;
class MGGelPositions;
class MGEdge;
class MGglViewAttrib;
class MGCommandDrawer;

#include <afxwin.h>
#include <gl\gl.h>
#include <list>
#include "mg/Types.h"
#include "mg/Pvector.h"
#include "mg/PickObjects.h"
#include "mg/Position.h"
#include "mg/Straight.h"
#include "mgGL/Color.h"
#include "mgGL/sysGLList.h"
#include "mgGL/ConstructionPlane.h"

/** @addtogroup DisplayHandling
 *  @{
 */

/// COpenGLWindow ビュー
class MGCLASS MGOpenGLView{

public:

///Compute eye position function.
typedef void (*Compute_eye_func)(
	void* calling_object,///<the class object instance pointer who uses this opglview.
	double eye_length,	///<the length from the center of the scene to eye position.
	MGPosition& eye,	///<Output eye position data computed from the above eye_length parameter.
	MGVector& view_up_vector///<view up vector when camera is located at eye.
);

//////////////////////////////PUBLIC MEMBER DATA/////////////

MGConstructionPlane m_cplane;///construction plane;

///System display list manager.
mgSysGLList m_sysgllist;

//////////////////////////////////METHOD/////////////////////

MGOpenGLView(
	bool perspective=true	///<indicates if the view is pespective or not.
);
///object, dvvf, and mmf must be set valid data before use of MGOpenGLView.

///Construct from MGglViewAttrib.
MGOpenGLView(
	const MGglViewAttrib& glatr
);

~MGOpenGLView();

/////////////////////////////////////////////////////////////////////////////////
public:

///Attach to MGCommandDrawer.
///Function's return value is the number of drawer's attached after attached.
size_t attach_drawer(MGCommandDrawer* drawer);

///Copy the informations of glview2 into this.
///Data that is not copied from glview2 are:
///m_hRC(Rendering context), m_hDC(Device Context), m_object, m_display_list,
///m_eye_func, and m_sysgllist.
///m_hRC and m_hDC can be set using setDCRC(), m_object can be set by set_object().
///m_eye_func by set_eye_func each.
///m_display_list can be made using make_display_list().
///m_sysgllist must be made by invoking openGL's
///display list generation functions.
void copy(const MGOpenGLView& glview2);

///copy the attributes in glatr into this.
void copy(const MGglViewAttrib& glatr);

///Get color values, Bcolor:background color, Gcolor:Object color,
///Hcolor:hilighted object color.
const MGColor& Bcolor()const{ return m_Bcolor;};
const MGColor& Gcolor()const{ return m_Gcolor;};
const MGColor& Hcolor()const{ return m_Hcolor;};

///return the center of the box of this view.
const MGPosition& center()const{return m_center;};
const MGPosition& center_current()const{return m_center_current;};

///clear current objects.
void clearCurrentObjects();

/// Clear all of the string displays.
///void clear_string_display();

const MGConstructionPlane& construction_plane()const{return m_cplane;};
MGConstructionPlane& construction_plane(){return m_cplane;};

///delete system display list by the function code.
bool DeleteDisplayList_by_function(size_t fc){
	return m_sysgllist.delete_lists_by_function_code(fc);
}

bool DeleteDisplayList_by_function_object_code(size_t fc,const MGGel* gel){
	return m_sysgllist.delete_lists_by_function_object_code(fc,gel);
}

///Detach and Return the detached MGCommandDrawer.
void detach_drawer(MGCommandDrawer* drawer);

///Obtainthe the diameter of the sphere that surround the whole model.
double diameter()const{return m_diameter;};

///Return display list name.
size_t display_list()const;

///Return line density for a surface to draw in wire mode.
int line_density()const{return m_line_density;};

///Draw the scene defined in this view including the current objects
///as hilighted.
///drawScene will invoke make_RC_current();
void drawScene(MGCL::VIEWMODE vmode, const MGPickObjects& pobjs);

/// grid snap.
void enable_grid_snap(bool bEnabled);

///Get the eye position.
const MGPosition& eye_position()const{return m_eyeP;}

///Free the openGL's context.
void FreeContext(){wglMakeCurrent( m_hDC, NULL ); }

///get ModelView matrix of OpenGL.
void get_model_matrix(
	double modelMat[16] ///<OpenGL's model matrix
)const;

///get projection matrix of OpenGL.
void get_projection_matrix(
	double projMat[16]	///<OpenGL's projection matrix
)const;

///Get the surface parameter value uv(u,v) where screen coordinate (sx,sy) is projected on.
///If no projection points are found, the nearest point to a perimeter of surf will
///be returned.
bool get_surface_parameter_glv(
	const MGFSurface& surf,
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).
	MGPosition& uv	///<surface parameter (u,v) where (sx,sy) is projected on will be returned.
)const;

///Get the current viewing environment.
///This output can be used for the parameters of initialize_viewing_environment,
///and thus the current environment that were save can be restored.
void get_viewing_environment(
	double& neard, double& fard,	///<near and far data
	double& maxsize, MGPosition& center,	///<Model's size(diameter) and the center
	double& scale, double& cx, double& cy,	///<Current viewing's scale, and panning data(cx, cy).
	MGPosition& eye, MGVector& up_vector,///<viewing eye position and view up vector
								///<These parameters are used fro gluLookAt. see the explanation.
	double currentMat[16]	///<OpenGL's ModelView current matrix.
)const;

///get viewport of OpenGL.
void get_viewport(
	int vp[4]	///<OpenGL's viewport will be output.
		///<(vp[0],vp[1]) is (left,bottom) coordinates.
		///<(vp[2],vp[3]) is (width, height) of the viewport.
)const;

///Get redering context.
HDC getHDC() const{return m_hDC;};

///Get redering context.
HGLRC getRC() const{return m_hRC;};

///Set the parent MGOpenGLView.
MGOpenGLView* get_parent_OpenGLView(){return m_parent_glView;};
const MGOpenGLView* get_parent_OpenGLView()const{return m_parent_glView;};
bool has_parent_OpenGLView()const{
	return m_parent_glView!=0;
}

///Test if this has a display list to draw.
bool has_display_list()const;

///Transform to the home position.
void home();

///Initialize the viewing environment.
///Initialization is done by the parameter box or all the necessary environmental data.
///initialize_viewing_environment will invoke make_RC_current(), and so, on return
///from initialize_viewing_environment, the rendering context is guaranteed to be current.
void initialize_viewing_environment(
	const MGBox& box,
	int view_num=1
		///<Standard view number:
		///<1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
		///<0: non standard view.
);
void initialize_viewing_environment(
	const MGglViewAttrib& viewA
);
void initialize_viewing_environment(
	double neard, double fard,	///<near and far data
	double maxsize, const MGPosition& center,	///<Model's size(diameter) and the center
	double scale, double cx, double cy,	///<Current viewing's scale, and panning data(cx, cy).
	const MGPosition& eye, const MGVector& up_vector,///<viewing eye position and view up vector
								///<These parameters are used fro gluLookAt. see the explanation.
	const double currentMat[16]	///<OpenGL's ModelView current matrix.
);
///Initialize the viewing environment to view objects projected onto a plane.
///This view is 2D view.
///center of the object is the origin of the plane.
void initialize_viewing_environment(
	const MGPlane& plane,
	double height,	///<height of the view,
					///<If rotate=true, height is the height after rotated.
	bool rotate=false///<if true, 90 degree rotation will be performed
					///<for the plane. (-v_deriv, u_deriv) will be (x,y),
					///<If false, (u_deriv, v_deriv) will be (x,y).
);

bool is_enabled_grid_snap()const;

///Return if this is a perspective view or not.
bool is_perspective() const{return m_perspective;};

///locate the screen coordinate (x,y) in the  3D world coordinate.
///(x, y)'s origin is (left, bottom) of the screen.
///In uv, the construction plane's parameter(u,v) coordinate will be returned.
MGPosition locate_glv(int x, int y, MGPosition* uv=0)const;

///construct the construction plane along with its display list.
void make_construction_plane(
	const MGBox& bx,
	int view_num=1	///<Standard view number:
	///<1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	///<0: non standard view.
);

void make_construction_plane(
	const MGPosition& mid,	///<center of the construction plane.
	const MGVector& uderi,	///<u-axis vector of the construction plane.
	const MGVector& vderi,	///<v-axis vector of the construction plane.
	double uspan,			///<span length between the lines along u-axis.
	double vspan,			///<span length between the lines along v-axis.
	size_t ulnum,			///<number of lines to draw along u-axis.
	size_t lnum			///<number of lines to draw along v-axis.
);

///Make openGL display list of a standard mgl file in the input glview.
///A standard mgl file means that all the objects are read into a MGGroup object,
///and that the objects are the target to make display list.
///If this is not the case, make display list as follows:
///1. Name the display list as e.g. NAME, and generate the display list by glNewName.
///2. Store all the data to display in the display list NAME.
///(Sometimes glCallList may be useful.)
///At this time, use glPushName() whose list name is size_t(MGObject*).
///MGOpenGLView's pick() will return the MGObject* when pick() is invoked.
///3. Let MGOpenGLView know the list NAME by invoking set_display_list().
size_t make_display_list(const MGGroup& grp);

/// Each viewport uses its own context, so we need to make sure the correct
/// context is set whenever we make an OpenGL command.
void make_RC_current()const{wglMakeCurrent( m_hDC, m_hRC ); }

///Translate and scale the current view.
///(x0, y0) to (x1,y1) is the rectangle of screen coordinate whose origin is
///(left,bottom).
void pan_zoom(int x0, int y0, int x1, int y1);

///Translate and scale the current view.
///box is the world coordinate's box cube.
void pan_zoom(const MGBox& box);

///Pick objects in the display list generated by make_display_list.
///Function's return value is MGPickObject vector.
///All the objects which were inside the pick aperture will be output.
///This data can be accessed using current_object()cy or current_PickObject().
///pick will invoke make_RC_current();
MGPickObjects pick_glv(
	MGCL::VIEWMODE vmode,///<View mode of the document.
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).

	double aperturex=-1.,///<specifies pick aperture of x and y.
	double aperturey=-1.,///<When <=0. value is specified, default value(the value
			///<obtained by pick_aperture() will be used.
	const MGAbstractGels& objtype=mgAll_Object,
			///<Target object kind. See MGGEL_KIND in "mg/types.h" or "mg/default.h"
	bool get_param=true
);

///Determine if screen coordinate (sx,sy) is closer to the start point or to the end
///of the curve curve.
///Functin's return value is 0: if start point, 1: if end point.
int pick_start_end_glv(
	const MGCurve& curve,
	int sx, int sy	///<Screen coordinates. (left, bottom) is (0,0).
);

///Pick a perimeter of the surface surf. That is, obtain the perimeter number
///that passes input (sx,sy) when drawn in the current view matrix.
///Function's return value is perimeter number picked.
///When no perimeters are picked, -1 will be returned.
int pick_perimeter_glv(
	const MGSurface& surf,
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv=0,	///<surface parameter (u,v) nearest to (sx,sy) will be returned.
	double aperturex=-1.,///<specifies pick aperture of x and y.
	double aperturey=-1.///<When <=0. value is specified, default value(the value
			///<obtained by pick_aperture() will be used.
);

///Pick an edge of the face f. That is, obtain the edge number
///that passes input (sx,sy) when drawn in the current view matrix.
///Function's return value is the edge pointer picked.
///When no edges are picked, null will be returned.
const MGEdge* pick_edge_glv(
	const MGFace& f,
	int sx, int sy,	///<Screen coordinates. (left, bottom) is (0,0).
	MGPosition* uv=0,	///<surface parameter (u,v) nearest to (sx,sy) will be returned.
	double aperturex=-1.,///<specifies pick aperture of x and y.
	double aperturey=-1.///<When <=0. value is specified, default value(the value
			///<obtained by pick_aperture() will be used.
);

///Get the pick aperture.
double pick_aperture()const{return m_pick_aperture;};

///Function's return value is the number of hit objects.
int MGOpenGLView::pick_to_select_buf(
	int sx, int sy,			///<Screen coordinates. (left, bottom) is (0,0).
	double aperturex, double aperturey,///<specifies pick aperture of x and y.
	size_t display_list,	///<display list that includes pick objects.
	size_t buf_size, GLuint* selectBuf,
		///<selected objects will be output in selectBuf inthe form of OpenGL.
	GLint* viewport=0,		///<view port will be output for the pick operation.
	GLdouble* modelMatrix=0,///<model matrix will be output for the pick operation.
	GLdouble* projMatrix=0	///<projection matrix will be output for the pick operation.
);

///project world coordinates to OpenGL's screen coordinates.
///If modelMat, projMat, or vp is not input, project will ask OpenGL to get them.
///Generally, users of project are recommended to get modelMat, projlMat, or
///vp, and input them to project.
///Before use of project, SetContext() must be invoked.
void project(
	const MGPosition& world, MGPosition& screen,
	const double* modelMat=0,	///<OpenGL's model matrix
	const double* projlMat=0,	///<OpenGL's projection matrix
	const int* vp=0				///<OpenGL's viewport
) const;

///Push back (function_code, object id) to the end or the beginning of
///system display list.
///Function's return value is the display list name to use for glNewList.
///sysgl can be handled(can be erased after the generation) by
/// (1)function code fc (2)object id oi(in other words, object name).
///If the system display list is to handle by the id(object name), oi must be
///meaningful. Usually oi is recommended to be the object pointer.
///If display list handling is not intended by the object id, set oi as null.
///size_t push_back_to_sysgl(size_t fc, MGGel* oi=0){return m_sysgllist.push_back(fc,oi);};
///size_t push_front_to_sysgl(size_t fc, MGGel* oi=0){return m_sysgllist.push_front(fc,oi);};
size_t push_back_to_sysgl(size_t fc,const MGGel* oi=0){return m_sysgllist.push_back(fc,oi);};
size_t push_front_to_sysgl(size_t fc,const MGGel* oi=0){return m_sysgllist.push_front(fc,oi);};

///Push back (function_code, object id) to the end or the beginning of
///system display list. sysgl must be a newed object, and the ownership will be
///transfered to this.
size_t push_back_to_sysgl(mgSysGL* sysgl){return m_sysgllist.push_back(sysgl);};

///Rotate the current view by the angle along the vector(x,y,z),
///performs a counterclockwise rotation of angle angle about
///the vector from the origin through the point (x, y, z).
void rotate(double angle, double x, double y, double z);

///Scale the current view by the factor.
void scale(double factor){m_scale*=factor;};

void set_center(const MGPosition& pos){ m_center = pos;}
void set_center_current(int x, int y);
void restore_center_current(){m_center_current=m_center;};

///Set display list MGGel*
void set_display_list(size_t list){m_display_list=list;};

///Set line density for a surface to draw in wire mode.
void set_line_density(int line_density=1){m_line_density=line_density;};

///Set the eye position computing function.
///this function is used wnly when the following function is invoked.
///"void initialize_viewing_environment(const MGBox& box);"
void set_eye_func(Compute_eye_func eye_func){m_eye_func=eye_func;};

///Set Projection matrix of OpenGL.
///set_frustum() will NOT invoke glLoadIdentity() first, so
///before use of set_frustum(), glLoadIdentity() must be invoked.
void set_frustum()const;

///Set Model_View matrix of OpenGL.
///set_model_matrix() will invoke glLoadIdentity() first.
void set_model_matrix();

///Set Object pointer who uses this MGOpenGLView.
///This pointer will be passed to Compute_eye_func
///as the 1st parameter as (*Compute_eye_func)(m_object,eye,m_eyeP).
void set_object(void* object){m_object=object;};

///Set the parent MGOpenGLView.
void set_parent_OpenGLView(MGOpenGLView* parent=0);

///Set if this view is a perspective view(true), or orthnormal view(falsle)
///, m_perspective.
void set_perspective(bool pers){m_perspective=pers;};

///Set the view-up vector.
void set_view_up_vector(const MGVector& up){m_up_vector=up;};

///Set background color;
void setBcolor(const MGColor& color);

///Set redering context.
///setDCRC will invoke make_RC_current();
void setDCRC(HDC dc, HGLRC rc);

void set_fovy(double fovy){m_fovy=fovy;};

void set_defalut_colors();

///Set default object color;
void setGcolor(const MGColor& color);

///Set hilight color;
void setHcolor(const MGColor& color);

void set_pick_aperture(double pick_aperture){m_pick_aperture=pick_aperture;};

///set the smooth factor of this view.
void set_smooth(double smooth){m_smooth=smooth;};

///Get the smooth factor of this view.
double smooth()const{return m_smooth;};

///Get the draw span length(approximate line segment length to draw curves).
double span_length()const{return m_diameter*smooth()*.3;};

///Display string at the position pos with the color
void draw_string(
	const CString& str,	///<String to display
	const MGPosition& pos,///<position for the string to display at.
	const MGColor* color=0///<color. When color=null, white color will be employed.
);

///SwapGLBuffers does not invoke make_RC_current() which is necessary
///before use of SwapGLBuffers.
void SwapGLBuffers(){::SwapBuffers(m_hDC);}

///Translate the current view by (dx, dy).
void translate(double dx, double dy);

///Translate the current view by (dx, dy) without current scale.
void translate_without_scale(double dx, double dy);

///Convert the windows screen coordinate (x,y) to MGCL's straight line.
///and get the intersection of the straight line and the construction plane.
///The origin of the screen coordinate is left, bottom. Not left, top.
///The direction of the sl is from the screen to the viewer.
void unproject(
	int x, int y,	///<screen coordinate whose origin is (left, bottom).
	MGStraight& sl,	///<The straight line of (x,y) will be returnred.
	MGCSisect& is	///<the intersectio of the sl and the construction plane
					///<will be returned.
)const;

///Convert the windows screen coordinate (x,y) to MGCL's straight line.
///The origin of the screen coordinate is left, bottom. Not left, top.
///The direction of the sl is from the screen to the viewer.
void unproject_to_sl_glv(int x, int y, MGStraight& sl)const;

///Get the view up vector.
const MGUnit_vector& view_up_vector()const{return m_up_vector;};

///Invoke wglUseFontBitmaps(getHDC(),0,STRING_COUNT,m_sitring_list_base);
bool UseFontBitmaps();

///compute the view volume far.
double view_volume_far()const{return m_far;}

///compute the view volume height.
double view_volume_height()const{return m_diameter/m_scale;}

///compute the view volume near.
double view_volume_near()const{return m_near;}

protected:
/// アトリビュート

///0.
	MGOpenGLView* m_parent_glView;///<If this OpenGLView has the parent, the pointer
		///<will be set. The parent means all of the display list are shared
		///<with the parent's, and on this drawScene's invocation, the parent's
		///<display list will be drawn.
	
	///Display list name of the document which is displayed on this glview.
	///This is used for glNewList display list name.
	///m_display_list=0 indicates no display list is generated for this glview.
	size_t m_display_list;

///1. Static atributes(usually not changed after this is initaialized)
	int m_line_density;///<line density for a surface to draw in wire mode.
	bool m_perspective;	///<Indicate if this is perspective or orthographic view.
						///<true if perspective.
	double m_fovy;		///<angle of top and bottom viewing pyramid in degrees.
						///<See gluPerspective's fovy.
						///<If m_perspective=false, m_fovy is not used.

	double m_near, m_far;///<Viewing frustum's near, far, and
	MGPosition m_eyeP;	///<eye data of gluLookAt.
	MGUnit_vector m_up_vector;///<up vector data of gluLookAt.

	///m_center and m_diameter define a sphere whose center and diameter are m_center and
	///m_diameter. The sphere includes the whole model.
	MGPosition m_center;///<World coordinate of the center of the document. 
	double m_diameter;	///<diameter of the sphere that sorround the model.

    HGLRC m_hRC;	///<Rendering context
	HDC m_hDC;		///<Device Context 

	MGColor m_Bcolor;	///<Background color.
	MGColor m_Gcolor;	///<Object lines color.
	MGColor m_Hcolor;	///<Object highlight color.

	double m_smooth;	///<Smoothness of the curves to draw.
	///< 1/smooth is the division number of a curve whose length is the window width.
	///< When smooth becomes small, smoothness increases.	
	double m_pick_aperture;///<Pick aperture. Number of pixels to allow picking.

	size_t m_sitring_list_base;///<wglUseFontBitmaps's base list number to display strings.

///2. Dynamic atributes(depend what is current viewing environment).
	MGPosition m_center_current;///<Current center after transformations.
	double m_scale;		///<Current scaling factor.

	///center of the screen in world coordinate
	///(when the center of the screen is supposed to be (0.,0.)).
	///Or (-m_cx, -m_cy) is the current panning translation distance.
	double m_cx, m_cy;

///3. Others.
	Compute_eye_func m_eye_func;///<function to compute eye position.
		///<this function is used wnly when the following function is invoked.
		///<"void initialize_viewing_environment(const MGBox& box);"

	void* m_object;	///<Generally this is a pointer who uses this MGOpenGLView.
		///<However, this is only used for Compute_eye_func
		///<to receive as the 1st parameter of their argumets:
		///< as (*Compute_eye_func)(m_object,eye,m_eyeP).
		
///Draw the scene defined in this view in wire mode, including the current objects as hilighted.
void drawWire(MGCL::VIEWMODE vmode,const MGPickObjects& pobjs);
		
///Draw the scene defined in this view in shading mode.
void drawShading();

private:

///Current command's picture drawer.
///MGOpenGLView will invoke the drawer last in the all the darwings.
std::list<MGCommandDrawer*> m_command_drawers;

///Convert the screen coordinate (sx, sy) to world coordinate (wx, wy) on the 
///view plane.
void screen_to_world(
	int wh[2],	///width(wh[0]) and height(wh[1]) of the screen.
	double sx,double sy, double& wx, double& wy
)const;

friend class COpenGLWindow;
friend class MGglViewAttrib;

};

///Convert the windows screen coordinate (x,y) to MGCL's straight line.
///The origin of the screen coordinate is left, bottom. Not left, top.
void unproject_to_sl(
	int x, int y,
	const double modelMatrix[16],	///<OpneGL's model matrix.
	const double projMatrix[16],	///<OpneGL's projection matrix.
	const int viewport[4],			///<OpneGL's viewport.
	MGStraight& sl
);

void get_near_position(
	const MGCurve* crv, 
	int sx,int sy,	///<screen coordinates whose origin is (left, bottom).
	const double modelMatrix[16],	///<OpneGL's model matrix.
	const double projMatrix[16],	///<OpneGL's projection matrix.
	const int viewport[4],			///<OpneGL's viewport.
	double& t	///parameter value of the curve crv near to (sx,sy) will be returned.
);

///////////////////////////////////////////////////////////////////////////////////////

/** @} */ // end of DisplayHandling group
#endif
