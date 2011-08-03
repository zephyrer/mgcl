/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGShell_HH_
#define _MGShell_HH_

#include <vector>
#include "mg/Default.h"
#include "mg/Pvector.h"
#include "mg/Position_list.h"
#include "mg/FPoint.h"
#include "topo/Complex.h"
#include "topo/Boundary.h"
#include "topo/CFisect_vector.h"
#include "topo/Face.h"
#include "topo/HHisect.h"
#include "topo/HHisect_vector.h"

class MGFSurface;
class MGIgesOfstream;

/** @addtogroup TOPO
 *  @{
 */

///MGShell is a composition of MGFace's(trimmed surface).
///See topology structure.
class MGCLASS MGShell: public MGBoundary{

public:

/////////Constructor/////////

///Default constructor.
MGShell(){;};

///Construct a shell of one face.
///Copy version of face.
MGShell(const MGFace& face);

///Construct a shell of one face.
///face is a pointer to newed face, and the ownership
///will be transfered to MGShell.
explicit MGShell(MGFace* face);

///Construct a shell of one face.
///Copy version of face.
MGShell(const MGFSurface& face);

///Construct a shell of one face.
///face is a pointer to newed face, and the ownership
///will be transfered to MGShell.
explicit MGShell(MGFSurface* face);

///Fundamental constructor.
///Construct from boundary complex(i.e. MGLoop).
///This constructor takes the ownership of MGCell* in boundary.
MGShell(
	std::list<MGCellNB*> boundaries
			///<Boundary data of the super class MGBoundary. List of faces.
);

/////////Destructor/////////

/////////operator overload/////////

///Assignment.
///When the leaf object of this and bnd2 are not equal, this assignment
///does nothing.
MGShell& operator=(const MGGel& gel2);
MGShell& operator=(const MGShell& gel2);

///Object transformation.
MGShell& operator+=(const MGVector& v){MGBoundary::operator+=(v);return *this;};
MGShell& operator-=(const MGVector& v){MGBoundary::operator-=(v);return *this;};
MGShell& operator*=(double scale){MGBoundary::operator*=(scale);return *this;};
MGShell& operator*=(const MGMatrix& mat){MGBoundary::operator*=(mat);return *this;};
MGShell& operator*=(const MGTransf& tr){MGBoundary::operator*=(tr);return *this;};

///comparison
bool operator<(const MGShell& gel2)const;
bool operator<(const MGGel& gel2)const;

///Debug Function

std::ostream& out(std::ostream& ostrm) const;

///IGES output function
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/////////Member Function/////////

///Test if this is an active boundary.
bool active() const{return false;};

///Append a face this is independent or connected to the already member faces of this
///shell. append_face() does not make any face join process.
///f must be a newed object, and the ownership will be transfered to this shell.
void append_face(MGFace* f){append_pcell(f);};

///Make a clone.
///Returned is pointer of newed object, must be deleted.
///When parent is specified, clone's parent is set to the parent.
MGShell* clone(MGCell& parent) const;
MGShell* clone() const;

///Make a clone that has not binders.
MGShell* clone_without_binders(MGCell& parent) const;
MGShell* clone_without_binders() const;

///Test if this is closed boundary.
bool closed() const{return false;};

///Compute the closest point from a point to this shell.
MGFPoint closest(const MGPosition& point) const;

///Delete a display list of this gel.
void delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	///<mgSysGL pointer will be apppended.
)const;

///////display member function.
virtual void display_arrows()const;
virtual void display_control_polygon()const;

//Ensure this shell has binder vertices of all the model edges.
//The binder vertices are stored in bcells of this shell(in MGComplex).
void ensure_BVertices_of_ModelEdges()const;

///Get the face pointer from its iterator in MGComplex of MGBoudarynD.
MGFace* face(pcellItr i);
const MGFace* face(const_pcellItr i)const;

///Get the face pointer from its pcell id in MGComplex of MGBoudarynD.
MGFace* face(size_t i);
const MGFace* face(size_t i)const;

/// Return This object's typeID
long identify_type() const;

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects intersection(const MGObject& obj2)const;
MGisects intersection(const MGCurve& obj2)const;
MGisects intersection(const MGFSurface& obj2)const;
MGisects intersection(const MGSurface& obj2)const;
MGisects intersection(const MGFace& obj2)const;
MGisects intersection(const MGShell& obj2)const;

///Intersection of a shell and a curve.
MGCFisect_vector isect(const MGCurve& curve) const;

///Intersection of two shells.
MGHHisect_vector isect(const MGShell& shell2)const;

///Intersection of a shell and a surface.
///This shell's face is face1 in HHisect and face2 is null.
MGHHisect_vector isect(const MGSurface& surf) const;

///Intersection of a shell and a face.
///This shell's face is face1 in HHisect and face2 is face.
MGHHisect_vector isect(const MGFace& face) const;

///Intersection of a shell and a face.
///This shell's face is face1 in HHisect and face2 is face.
MGHHisect_vector isect(const MGFSurface& face) const;

///Make 2 types of display list of this gel(wire and shading).
///Return is the display list name.
size_t make_display_list(
	double span_length,///span length to approximate by polyline.
	int line_density=1///line density to draw surface in wire mode.
)const;

///Make a display list without color of this gel.
///Return is the display list name.
size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Get manifold dimension.
unsigned manifold_dimension() const{return 2;};

///Merge a face at free edges of this shell.
///Function's return value is
///   false: merge was not done because no common edges were found.
///   true: merge of the shell and the face was done.
///      Case of true includes that merge was done only partialy.
///      This case occurs when after some common edges are merged,
///      some common edges are found to have contradictionary direction.
///      Merge is done only for the first common edges found to have
///      same direction.
///When the function's return value is false, the face will not be added to
///the shell. Merge negates the input face direction when this shell's direction
///is not the same as the input face's one along the common edge(s).
///The second form is a pointer version. The face must be newed object
///and will be destructed after merge when function's return value is true
///(merge was processed).
///This shell may be dummy shell. In this case, the face is added to the 1st
///face in the shell.
bool merge_at_common_edge(const MGFace& face);
bool merge_at_common_edge(MGFace* face);
bool merge_at_common_edge(const MGFSurface& face);
bool merge_at_common_edge(MGFSurface* face);

///Get the number of faces included in this shell.
size_t number_of_faces()const{return number_of_pcells();};

///Test if a point is on the shell or not.
///Function's return value is true if the point is on the shell, and false if not.
///The point parameter of the shell is returned in fp if true is returned.
///If false is returned, the closest point of the shell will be returned in fp.
bool on(const MGPosition& point,
	MGFPoint& fp				///<Shell's point parameter value.
	) const;

///Obtain perpendicular points of a shell from a point.
std::vector<MGFPoint> perps(const MGPosition& point) const;

///Obtain the projected curves of a curve onto the shell.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the shell if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves(iline() of the MGHHisect members of lines),
///and the other is (u,v) curves of the parameter space of the surfaces
///(uvline1() of the MGHHisect members of lines ).
///*** uvline2() of the MGHHisect members of lines is a deque of length zero
///(not used).
///Function's return value is the number of curves obtained.
///When <0 is returned, some internal error occured.
int project(
	const MGCurve& crv,		///<given world coordinate curve to project.
	MGHHisect_vector& lines,
			///<World coordinates (x,y,z) lines and (u,v) lines of the projected
			///<curves will be returned in lines.
	const MGVector& vec=mgNULL_VEC	///<projection direction.
			///<if vec = NULL, then projection that is normal to the shell.
	)const;

///Shade the object in world coordinates.
void shade(
	double span_length	///<Line segment span length.
)const;

///Return MGShell pointer if this MGGel is an MGShell, else return null.
MGShell* shell(){return this;};
const MGShell* shell()const{return this;};

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
MGPvector<MGCurve> skeleton(int density=1)const;

///Obtain all the parameter curves at knots of u and v knot vector.
MGPvector<MGCurve> skeleton_at_knots()const;

virtual std::string whoami()const{return "Shell";};

protected:

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:

///Copy constructor with mapping.
///Binder cells of the pcells in loop will be registered in cmap
MGShell(
	const MGShell& shell,	///<original shell.
	MGCellMap& cmap);		///<cellmap to register binder association.

///Make a clone.
///The forms that have cmap as an argumetnt is to register binder association.
///Binder cells of the pcells in this boundary will be registered in cmap.
///Returned is pointer of newed object, must be deleted.
///When parent is specified, clone's parent is set to the parent.
MGShell* clone(MGCell& parent, MGCellMap& cmap) const;
MGShell* clone(MGCellMap& cmap) const;

///Get the partner point uvuv_id of the uvuv_list[j].
///If found, return true, if not, return false.
///Shell*(face or surface) version.
///The other one must not be shell, a surface or a face.
bool isect_partner(
	const MGPosition& uvuv,///<The point to get the partner.
	MGPosition_list* uvuv_list,
	size_t& j,		///<Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr)const;
		///<Partner point's iterator of uvuv_list[j] will be returned.

///Get the partner point uvuv_id of the uvuv_list[j].
///If found, return true, if not, return false.
///Shell*shell intersection version.
bool isect_partner(
	const MGPosition& uvuv,///<The point to get the partner.
	MGPosition_list* uvuv_list,
	size_t& i,
	size_t& j,		///<Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr)const;
		///<Partner point's iterator of uvuv_list[i+nf1*j] will be returned.

///Make a display list of only call lists of this shell's faces.
void make_only_call_list(
	int mode	///<  =0 when wire, =1 when wire and no color, =2 when SHADING.
)const;

};

/** @} */ // end of TOPO group
#endif
