/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGel_HH_
#define _MGGel_HH_

#include <iosfwd>
#include <string>
#include "mg/MGCL.h"
#include "mg/AbstractGels.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"
#include "mg/Pvector.h"
#include "mg/Transf.h"

class MGAttribedGel;
class MGAttrib;
class MGIfstream;
class MGOfstream;
class MGAttrib;
class MGGroup;
class MGObject;
class MGGeometry;
class MGPoint;
class MGCurve;
class MGSurface;
class MGTopology;
class MGFace;
class MGShell;
class MGOpenGLView;
class MGIgesOfstream;
class mgSysGL;

/** @defgroup GelRelated Gel Related class
 *  MGGel is top abstract class for MGObject, MGGroup, and MGGLAttrib.
 *  @{
 */

//Define MGGel Class.

///MGGel is an abstract class which represents a group element.
///Gel is the abbreviation of group element, is designed to store
///in MGGroup as an element.
///Subclasses of MGGel are:
///(1) MGAttribedGel(whose sub are MGObject, MGGroup), or (2) MGAttrib.
///MGGel provides functions of serialization of objects.
///All the objects of MGGel subclasses can be serialized using
///MGGroup::make_file(), and MGGroup constructor.
class MGCLASS MGGel{

///string stream function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGGel& gel);

public:

////////////Constructor////////////

//Void constructor(初期化なしでオブジェクトを作成する。)
//MGGel();

//Copy constructor.
//MGGel(const MGGel& obj2);

///Virtual Destructor
virtual ~MGGel();

///Assignment.
///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
virtual MGGel& operator=(const MGGel& gel2){return *this;};

///Comparison.
virtual bool operator==(const MGGel& gel2)const{return false;};
virtual bool operator!=(const MGGel& gel2)const{return !(operator==(gel2));};
virtual bool operator<(const MGGel& gel2)const;
virtual bool operator>(const MGGel& gel2)const{return gel2<(*this);};

////////////Member Function////////////

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const=0;

///IGES output function
///(Default function is no operation to output)
///Function's return value is the directory entry id created.
virtual int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const{return 0;};

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
virtual MGGel* clone()const=0;

///Delete a display list of this gel.
virtual void delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	///<mgSysGL pointer will be apppended.
)const{;};

///Obtain display list name. 0(null) means this gel need not to be displayed.
///The default is the pointer of this gel.
virtual size_t dlist_name()const{return size_t(this);};

///draw attribute data.
virtual void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///Make a display list of this gel.
///Return is the display list name.
virtual size_t make_display_list(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const{return 0;};

///Make a display list without color of this gel.
///Return is the display list name.
virtual size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const{return 0;};

///Return MGAttrib pointer if this MGGel is an MGAttrib, else return null.
virtual MGAttrib* attrib(){return 0;};
virtual const MGAttrib* attrib()const{return 0;};

///Return MGGroup pointer if this MGGel is an MGGroup, else return null.
virtual MGGroup* group(){return 0;};
virtual const MGGroup* group()const{return 0;};

///Return MGObject pointer if this MGGel is an MGObject, else return null.
virtual MGObject* object(){return 0;};
virtual const MGObject* object()const{return 0;};

///Return MGGeometry pointer if this MGGel is an MGGeometry, else return null.
virtual MGGeometry* geometry(){return 0;};
virtual const MGGeometry* geometry()const{return 0;};

///Return point pointer if this MGGel is an MGPoint, else return null.
virtual MGPoint* point(){return 0;};
virtual const MGPoint* point()const{return 0;};

///Return curve pointer if this MGGel is an MGCurve, else return null.
virtual MGCurve* curve(){return 0;};
virtual const MGCurve* curve()const{return 0;};

///Return MGSurface pointer if this MGGel is an MGSurface, else return null.
virtual MGSurface* surf(){return 0;};
virtual const MGSurface* surf()const{return 0;};

///Return MGTopology pointer if this MGGel is an MGTopology, else return null.
virtual MGTopology* topology(){return 0;};
virtual const MGTopology* topology()const{return 0;};

///Return MGFace pointer if this MGGel is an MGFace, else return null.
virtual MGFace* face(){return 0;};
virtual const MGFace* face()const{return 0;};

///Return MGShell pointer if this MGGel is an MGShell, else return null.
virtual MGShell* shell(){return 0;};
virtual const MGShell* shell()const{return 0;};

/// Return This object's typeID
virtual long identify_type() const = 0;

///Test if this gel includes an object.
virtual const MGObject* includes_object()const=0;
virtual MGObject* includes_object()=0;

///Test if this gel should be displayed or not.
///True: not to display, false:to display.
virtual bool no_display()const{return false;};

///Output the content as std::string.
///The output string is the same as cout<<MGGel.
std::string string_content()const;

///Transform the gel by the argument.

///translation
virtual void transform(const MGVector& v){;};

///scaling.
virtual void transform(double scale){;};

///matrix transformation.
virtual void transform(const MGMatrix& mat){;};

///general transformation.
virtual void transform(const MGTransf& tr){;};

///Determine if this is one of the input types or not.
///Function's return value is true if this is one of the input types.
bool type_is(const MGAbstractGels& types)const;

///////display member function.
virtual void display_arrows()const{;};
virtual void display_break_points()const{;};
virtual void display_control_polygon()const{;};
virtual void display_curvatures(
	double	scale,	///<scaling of the graph.
	int		density,///<densitiy of the graph.
	bool	use_radius///<true:radius display, false:curvature display.
)const{;};

virtual std::string whoami()const=0;

protected:

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

private:

friend class MGIfstream;
friend class MGOfstream;

};

///Construct a null newed MGGel from the type id TID.
///Objects handled by MGIfstream or MGOfstream is only the following objects.
MGDECL MGGel* MGNullGel(long TID);

/** @} */ // end of GelRelated group
#endif
