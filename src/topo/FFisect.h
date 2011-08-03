/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFFisect_HH_
#define _MGFFisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Curve.h"
#include "mg/FPline.h"
#include "topo/Face.h"

// MGFFisect.h
// Header for MGFFisect

//Forward Declaration
class MGObject;
class MGCurve;
class MGFace;

///MGFFisect represents one intersection line of a MGFace and MGFace or MGSurface..
///The behavior of MGFFisect is like a auto_ptr. Copy or assignment
///of MGFFisect means transfer of the ownership of all the included curve
///to copied or assigned MGFFisect and original MGFFisect does not have the
///ownership of the curves any more. User should be aware of it.
class MGCLASS MGFFisect:public MGisect{

public:

//////////// Constructor ////////////

///Void constructou. èâä˙âªÇ»ÇµÇ≈Dummyåê¸Çê∂ê¨
MGFFisect():m_curve(0){;};

///Construct providing all the raw data.
///The ownership of iline, face1, and face2 are all transfered to MGFFisect.
///All of these objects must be newed ones.
MGFFisect(
	MGCurve* iline,	///<Pointer of a newed object.
	const MGFace* face1,
	MGCurve* param1,///<Pointer of a newed object.
	const MGFace* face2,
	MGCurve* param2///<Pointer of a newed object.
):m_curve(iline),m_face1line(face1,param1),m_face2line(face2,param2){;};

///Construct providing all the raw data.
///Copy version. Copy of the three curves will take place.
MGFFisect(
	const MGCurve& iline,
	const MGFace* face1,
	const MGCurve& param1,
	const MGFace* face2,
	const MGCurve& param2
);

///Construct providing all the raw data.
MGFFisect(
	MGCurve* iline,	///<Pointer of a newed object.
	MGFPline face1uv,
	MGFPline face2uv
);

/// Copy Constructor;
/// ffi's ownership of all the three curves will be released.
MGFFisect(const MGFFisect& ffi);

//////////// Destructor ////////////
///~MGFFisect();

//////////// Operator overload ////////////

///Assignment
/// ssi's ownership of all the three curves will be released.
MGFFisect& operator= (const MGFFisect& ffi);

///Comparison operator.
bool operator< (const MGFFisect& ssi2)const;
bool operator> (const MGFFisect& ssi2)const{return ssi2<(*this);};
bool operator<= (const MGFFisect& ssi2)const{return !(ssi2<(*this));};
bool operator>= (const MGFFisect& ssi2)const{return !((*this)<ssi2);};
bool operator== (const MGFFisect& ssi2)const;
bool operator!= (const MGFFisect& ssi2)const{return !operator==(ssi2);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
bool operator< (const MGCFisect& is)const{return false;};
bool operator< (const MGSSisect& is)const{return false;};
bool operator== (const MGisect& is)const;

//////////// Memeber Function ////////////

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return *m_curve;};

///Return the 1st object's parameter value curve of the intersection.
const MGCurve* isect1_param1()const{return &(m_face1line.uvline());};

///Return the 2nd object's parameter value curve of the intersection.
const MGCurve* isect1_param2()const{return &(m_face2line.uvline());};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 1;};

/// Output function.
std::ostream& out(std::ostream& ostrm)const;

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12();

protected:

	///Get the object1 pointer.
	const MGObject* object1(const MGObject* obj)const{
		return m_face1line.face()->object_pointer();
	};

	///Get the object2 pointer.
	const MGObject* object2(const MGObject* obj)const{
		return m_face2line.face()->object_pointer();
	};

private:

//////////// Member Data ////////////
	
	MGCurve* m_curve;	///<world coordinate curve expression. A newed object.
	MGFPline m_face1line;///<1st Face's parameter line expression.
	MGFPline m_face2line;///<2nd Face's parameter line expression,
						///<may be null when the objective is not a face but a surface.
};

/** @} */ // end of IsectContainer group
#endif
