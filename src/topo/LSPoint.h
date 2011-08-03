/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLSPoint_HH_
#define _MGLSPoint_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Vector.h"

//
//Define MGLSPoint Class.

class MGEdge;
class MGPosition;

///MGLSPoint is to express a loop and a surface intersection point.
///The expression is {MGEdge* binder, double tb,(u,v)}, where binder is
///binder edge of the loop, tb is parameter value of the binder, and (u,v) is
///the surface parameter value.
class MGCLASS MGLSPoint{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLSPoint& );

/////////Constructor/////////
MGLSPoint():m_pedge(0){;};

///Construct from all the necessary data.
MGLSPoint(
	const MGEdge* pedge,///<Loop's parameter edge pointer.
	double t,		///<Parameter values of pedge's binder edge's curve.
	double u, double v)///<Surface parameter value.
	:m_pedge(pedge), m_t(t),m_u(u), m_v(v){;};

/////////Operator oveload/////////

///Comparison operator.
bool operator< (const MGLSPoint& ls2)const;
bool operator> (const MGLSPoint& ls2)const{return ls2<(*this);};
bool operator<= (const MGLSPoint& ls2)const{return !(ls2<(*this));};
bool operator>= (const MGLSPoint& ls2)const{return !((*this)<ls2);};
bool operator== (const MGLSPoint& ls2)const;
bool operator!= (const MGLSPoint& ls2)const{return !operator==(ls2);};

/////////Member function/////////

///Return the binder edge pointer.
const MGEdge* parameter_edge()const{return m_pedge;};

///Return binder edge's parameter data.
double binder_param()const{return m_t;};

///Set binder edge pointer.
void set_pedge(const MGEdge* pedge){ m_pedge=pedge;};

///Set binder edge pointer.
void set_binder_param(double t){ m_t=t;};

///Set binder edge pointer.
void set_surface_param(double u, double v){m_u=u; m_v=v;};
void set_surface_param(const MGPosition& uv);

///Return surface's parameter data.
MGPosition surface_param()const;
void surface_param(double& u, double& v)const;

///Obtain world point coordinate data from the binder edge.
MGVector world_point()const;

private:
	const MGEdge* m_pedge;	///<Parameter edge pointer of the ip.
	double m_t;			///<m_pedge's curve parameter value.
		///<(m_t is not the parameter of m_pedge, but of m_pedge->binder_edge().
	double m_u, m_v;	///<Surface parameter value (u,v).

};

/** @} */ // end of IsectContainer group
#endif
