/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCSisect_HH_
#define _MGCSisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Point.h"

// MGCSisect.h
// Header for MGCSisect

//Forward declaration.
class MGInterval;
class MGPosition;

///
///One Intersection of curve and surface.
///If more than one are necessary to hold, CSisect_list should be used.
/// Curve と Surface の交点を表現する。
class MGCLASS MGCSisect:public MGisect{

public:

//////////////Constructor////////////

///Void Constructor(初期化なしで交点を生成)
MGCSisect ();

///Input all necessary components(全てのコンポーネントを指定して交点を生成)
MGCSisect(
	const MGPosition& point,	///<intersection point.
	double t,								///<Curve's parameter value.
    const MGPosition& uv,					///<Surface's parameter values.
	const MGCSRELATION rl=MGCSREL_UNKNOWN	///<Curve and Surface relation
);

////////////Operator overload////////////

bool operator< (const MGCSisect& csi)const;
bool operator> (const MGCSisect& csi)const{return csi<(*this);};
bool operator<= (const MGCSisect& csi)const{return !(csi<(*this));};
bool operator>= (const MGCSisect& csi)const{return !((*this)<csi);};
bool operator== (const MGCSisect& csi)const;
bool operator!= (const MGCSisect& csi)const{return !operator==(csi);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCFisect& is)const{return true;};
bool operator< (const MGSSisect& is)const{return true;};
bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

////////////Member Function////////////

///obtain the distance in parameter space.
void distance(const MGCSisect& isect2,	///2nd isect.
			  double& t, double& u, double& v)const;

///Compute square of parameter space distance between this and
///isect2.
double distance_square(const MGCSisect& isect2)const;

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12(){;};

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return m_ipoint;};

///Return the 1st object's parameter value of the intersection.
MGPosition isect0_param1()const{return MGPosition(1,&m_t);};

///Return the 2nd object's parameter value of the intersection.
MGPosition isect0_param2()const{return m_uv;};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 0;};

/// Output virtual function.
std::ostream& out(std::ostream& ostrm)const;

///Return parameter value of curve.
/// 交点の curve のパラメータ値を返却する。
double param_curve() const{return m_t;};

///Return parameter value of surface.
/// 交点の Surface のパラメータ値を返却する。
const MGPosition& param_surface() const{return m_uv;};

///Return coordinate values(交点の座標値を返却する)
const MGPosition& point() const{return m_ipoint.position();};

///Return Surface and curve relationship at the i.p.
///交点での関係を返却する。
MGCSRELATION rel() const{return m_rel;};

private:
///Member data
	
	MGPoint m_ipoint; ///< coordinate values(交点の座標値)
	double m_t;	         ///< parameter value of curve
						 ///<(交点におけるcurve のパラメータ値)
	MGPosition m_uv;     ///< parameter value of surface
						 ///<(交点におけるsurface のパラメータ値)
	MGCSRELATION m_rel;	 ///< Surface and curve relationship at the i.p.
						 ///< 交点におけるCurveとSurfaceとの関係

};

/** @} */ // end of IsectContainer group
#endif
