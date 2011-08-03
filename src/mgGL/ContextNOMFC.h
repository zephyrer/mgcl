/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGContext_HH_
#define _MGContext_HH_

#include "mg/Pvector.h"
#include "mg/Attrib.h"
#include "mgGL/Color.h"
#include "mgGL/Appearance.h"
//#include "mgGL/SnapAttrib.h"
//#include "mgGL/glViewAttrib.h"
#include "Tl/TLInputParam.h"

///////////////////////////////////////////////////

class MGOpenGLView;
class MGOfstream;
class MGIfstream;

/** @addtogroup DisplayHandling
 *  @{
 */

/// MGContext defines the attributes of a document.
class MGCLASS MGContext: public MGAttrib{

public:

MGContext():m_line_density(1),m_smooth(float(.01)),
m_pick_aperture(float(5.)),m_appearance(0),m_gel(0),
m_Bcolor(.5,.5,.5),//Background color.
m_Gcolor(0.,0.,0.),//Object lines color.
m_Hcolor(1.,1.,0.)//Object highlight color.
{	//for(int i=0; i<4; i++) m_views[i]=0;
	m_tolerance[0]=0.5E-3;			//wc_zero
	m_tolerance[1]=1.0E-6;			//rc_zero
	m_tolerance[2]=1.0E-20;			//mach_zero
	m_tolerance[3]=m_tolerance[0];	//line_zero;
	m_tolerance[4]=.0025;			//angle_zero;
	m_tolerance[5]=5.0E+2;			//max_knot_ratio;
};

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
MGContext* clone()const{return new MGContext(*this);};

/// set/get snap attrib data.
//const MGSnapAttrib& snap_attrib()const{return m_snap_attrib;}; //Snap data
//MGSnapAttrib& snap_attrib(){return m_snap_attrib;}; //Snap data
//void set_snap_attrib(const MGSnapAttrib& snap_attrib){m_snap_attrib=snap_attrib;};

//set/get line_density attrib data.
int line_density()const{return m_line_density;};
void set_line_density(int line_density){m_line_density=line_density;};

/// set/get color data.

///Background color.
const MGColor& Bcolor()const{return m_Bcolor;};

///Object lines color.
const MGColor& Gcolor()const{return m_Gcolor;};	

///Object highlight color.
const MGColor& Hcolor()const{return m_Hcolor;};	

///Background color.
MGColor& Bcolor(){return m_Bcolor;};	

///Object lines color.
MGColor& Gcolor(){return m_Gcolor;};	

///Object highlight color.
MGColor& Hcolor(){return m_Hcolor;};	

///set appearance.
void set_appearance(
	MGAppearance* appearance	///<appearance must be a newed object. The ownership will
								///<be transfered to this MGContext.
								){;};

/// set/get smooth data.
///The smooht data is used for MGOpenGLView. See the explanation.
///Smoothness of the curves to draw.
double smooth()const{return m_smooth;};	

void set_smooth(double smooth){m_smooth=float(smooth);};

///pick_aperture.
///The pick_aperture data is used for MGOpenGLView. See the explanation.
double pick_aperture()const{return m_pick_aperture;};
void set_pick_aperture(double pick_aperture){m_pick_aperture=float(pick_aperture);};

/// Tolerance data.
double* tolerance(){return m_tolerance;};
const double* tolerance()const{return m_tolerance;};

/// Tessellation parameter.
const mgTLInputParam& tessellate_param()const{return m_tessellate_param;};
mgTLInputParam& tessellate_param(){return m_tessellate_param;};
void set_tessellate_param(const mgTLInputParam& tessellate_param){
	m_tessellate_param=tessellate_param;};

/// Appearance data.
const MGAppearance* appearance()const{return m_appearance;};
MGAppearance* appearance(){return m_appearance;};

/// Return This object's typeID
long identify_type() const{return MGCONTEXT_TID;};

/////////Attributes execution functions.///////

///stream text output.
std::ostream& out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"MGContext="<<this<<std::endl;
	return ostrm;
};

std::string whoami()const{return "Context";};

///Read all member data.
void ReadMembers(MGIfstream& buf){;};

///Write all member data
void WriteMembers(MGOfstream& buf)const{;};

private:

	//MGSnapAttrib m_snap_attrib;///<Snap data
	
	///////Folowing are MGOpenGLView attributes./////////

		int m_line_density;///<line density for a surface to draw in wire mode.
		MGColor m_Bcolor;	///<Background color.
		MGColor m_Gcolor;	///<Object lines color.
		MGColor m_Hcolor;	///<Object highlight color.
		float m_smooth;	///<Smoothness of the curves to draw.
			///< 1/smooth is the division number of a curve whose length is the window width.
			///< When smooth becomes small, smoothness increases.	
		float m_pick_aperture;///<Pick aperture. Number of pixels to allow picking.
		//MGglViewAttrib* m_views[4];///<Newed objects of MGglViewAttrib, may be null.
			///<[0]=for 3D perspective view
			///<[1]=for (x,y) 2D view
			///<[2]=for (y,z) 2D view
			///<[3]=for (z,x) 2D view

	///////Above are MGOpenGLView attributes./////////

	double m_tolerance[6];///<MGCL Tolerance.
			///<[0]=wc_zero;
			///<[1]=rc_zero;
			///<[2]=mach_zero;
			///<[3]=line_zero;
			///<[4]=angle_zero;
			///<[5]=max_knot_ratio;
	mgTLInputParam m_tessellate_param;///<tessellation parameter.
	MGAppearance* m_appearance;///<Newed objects of MGAppearance.

	MGGel* m_gel;				///<Reserved area for future extension.
};

/** @} */ // end of DisplayHandling group
#endif // _MGContext_HH_
