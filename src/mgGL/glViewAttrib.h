/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGLViewAttrib_HH_
#define _MGGLViewAttrib_HH_

#include "mg/Position.h"
#include "mgGL/ConstructionPlane.h"
#include <iosfwd>

class MGOfstream;
class MGIfstream;
class MGBox;
class MGOpenGLView;

/** @addtogroup DisplayHandling
 *  @{
 */

///MGglViewAttrib is a class to serialize MGOpenGLView, used by MGContext.
class MGCLASS MGglViewAttrib{

public:
	
///Debug Function.
MGDECL friend std::ostream& operator<< (std::ostream& out, const MGglViewAttrib& atr);

/// Serialization fucntion.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGglViewAttrib& atr);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGglViewAttrib& atr);

//////////////////////////////////METHOD/////////////////////
MGglViewAttrib(){;};

MGglViewAttrib(
	const MGOpenGLView& glview
);

///~MGglViewAttrib();

/////////////////////////////////////////////////////////////////////////////////

///Copy the informations of glview2 into this.
void copy(const MGOpenGLView& glview);

///return the center of the box of this view.
const MGPosition& center()const{return m_center;};
MGPosition& center(){return m_center;};

const MGConstructionPlane& construction_plane()const{return m_cplane;};
MGConstructionPlane& construction_plane(){return m_cplane;};

///Obtainthe the diameter of the sphere that surround the whole model.
double diameter()const{return m_diameter;};
double& diameter(){return m_diameter;};

///Get the eye position.
const MGPosition& eye_position()const{return m_eyeP;}
MGPosition& eye_position(){return m_eyeP;}

///Return if this is a perspective view or not.
bool is_perspective() const{return m_perspective;};

///set the viewing environment info.
void set(
	double neard, double fard,	///<near and far data
	double diameter, const MGPosition& center,	///<Model's size(diameter) and the center
	double scale, double cx, double cy,	///<Current viewing's scale, and panning data(cx, cy).
	const MGPosition& eye, const MGVector& up_vector,///<viewing eye position and view up vector
								///<These parameters are used fro gluLookAt. see the explanation.
	const double currentMat[16]	///<OpenGL's ModelView current matrix.
);

///Set if this view is a perspective view(true), or orthnormal view(falsle)
///, m_perspective.
void set_perspective(bool pers){m_perspective=pers;};

///Set the view-up vector.
void set_view_up_vector(const MGVector& up){m_up_vector=up;};

///Get the view up vector.
const MGVector& view_up_vector()const{return m_up_vector;};
MGVector& view_up_vector(){return m_up_vector;};

///compute the view volume far.
double view_volume_far()const{return m_far;}

///compute the view volume height.
double view_volume_height()const{return m_diameter/m_scale;}

///compute the view volume near.
double view_volume_near()const{return m_near;}

private:
/// アトリビュート
	bool m_perspective;	///<Indicate if this is perspective or orthographic view.
						///<true if perspective.
	double m_fovy;		///<angle of top and bottom viewing pyramid in degrees.
						///<See gluPerspective's fovy.

	double m_near, m_far;///<Viewing frustum's near, far, and
	MGPosition m_eyeP;	///<eye data of gluLookAt.
	MGVector m_up_vector;///<up vector data of gluLookAt.

	///m_center and m_diameter define a sphere whose center and diameter are m_center and
	///m_diameter. The sphere includes the whole model.
	MGPosition m_center;///<World coordinate of the center of the document. 
	double m_diameter;	///<diameter of the sphere that sorround the model.

	double m_scale;		///<Current scaling factor.
	double m_cx, m_cy;	///<center of the screen in world coordinate
			///<(when the center of the screen is supposed to be (0.,0.).
			///<or (-m_cx, -m_cy) is the current panning translation distance.

	MGConstructionPlane m_cplane;///<construction plane;

	double m_currentMat[16];

friend class MGOpenGLView;

};

/** @} */ // end of DisplayHandling group
#endif
