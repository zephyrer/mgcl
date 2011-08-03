/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCFisect_HH_
#define _MGCFisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/CSisect.h"

class MGPosition;
class MGFace;

//
//Define MGCFisect Class.

///MGCFisect is to represent an intersection of a face and a curve.
///(MGCSisect csi, MGFace* f) where csi consists of world point, curve parameter,
///and face(surface) parameter, and f is a face pointer.
class MGCLASS MGCFisect:public MGisect{

public:
/////////Constructor/////////

///void constructor.
MGCFisect():m_face(0){;};

///Construct from all the necessary data.
MGCFisect(
	const MGCSisect& csi,	///<isect data (point, curve parameter value,
							///<            surface parameter value)
	const MGFace& face		///<face.
):m_csi(csi), m_face(&face){;};

///Construct from all the necessary data.
MGCFisect(
	const MGPosition& point,	///<World coordinate point data of the isect.
	const double& t,			///<curve parameter value of the isect.
	const MGPosition& uv,		///<Face(Surface) parameter value of the isect.
	const MGFace& face		///<face.
);

/////////Operator oveload/////////

bool operator< (const MGCFisect& fp)const;
bool operator> (const MGCFisect& fp)const{return fp<(*this);};
bool operator<= (const MGCFisect& fp)const{return !(fp<(*this));};
bool operator>= (const MGCFisect& fp)const{return !((*this)<fp);};
bool operator== (const MGCFisect& fp)const;
bool operator!= (const MGCFisect& fp)const{return !operator==(fp);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
bool operator< (const MGSSisect& is)const{return true;};
bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

/////////Member function/////////

///Return isect data.
const MGCSisect& csi()const{return m_csi;};

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12(){;};

///return the face.
const MGFace& face()const{return *m_face;};

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return m_csi.isect();};

///Return the 1st object's parameter value of the intersection.
MGPosition isect0_param1()const{return m_csi.isect0_param1();};

///Return the 2nd object's parameter value of the intersection.
MGPosition isect0_param2()const{return m_csi.isect0_param2();};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 0;};

/// Output virtual function.
std::ostream& out(std::ostream& ostrm)const;

///Return coordinate values(交点の座標値を返却する)
const MGPosition& point() const{return m_csi.point();};

///Return the parameter value of the curve.
/// 交点の curve のパラメータ値を返却する。
double param_curve() const{return m_csi.param_curve();};

///Return the parameter value of the surface.
/// 交点の Surface のパラメータ値を返却する。
const MGPosition& param_face() const{return m_csi.param_surface();};

///Return Surface and curve relationship at the i.p.
///交点での関係を返却する。
///MGCSRELATION rel() const{return m_csi.rel();};

private:
	MGCSisect m_csi;		///<curve and surface intersection, includes
							///<world coordinates, curve parameter, and surface parameter
							///<of the intersection point.
	const MGFace* m_face;	///<face pointer of the intersection point.

};

/** @} */ // end of IsectContainer group
#endif
