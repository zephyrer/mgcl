/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// CurveContinuity.h

#if !defined( __MGCurveContinuity_H__)
#define __MGCurveContinuity_H__

#include "mg/Unit_vector.h"
#include "mg/Position.h"
class MGCurve;

/** @addtogroup GEORelated
 *  @{
 */

/// Curve continuity measuring class.
/// MGCurveContinuity measures the continuity of two curves.
/// Measuring is done at the closest end points of the two curves.
///To us MGCurveContinuity, construct MGCurveContinuity object by inputting
///two curves. All of the continuity information will be generated in the
///MGCurveContinuity object.
class MGCurveContinuity{
public:
	/// enumeration to represent geometric continuity
	enum CONTINUITY{
		DISCONT=-1,	///< means not continuous
		G0=0,		///< means G0 continuity
		G1,			///< means G1 continuity
		G2,			///< means G2 continuity
	};

	///constructor
	MGCurveContinuity(
		const MGCurve& curve1,
		const MGCurve& curve2
	);

	///Get the continuity of the two curves.
	CONTINUITY get_continuity()const{return m_continuity;};

	///Get the curve1 position(start or end) that is closest to curve2.
	const MGPosition& P1()const{return m_P1;};

	///Get the curve2 position(start or end) that is closest to curve1.
	const MGPosition& P2()const{return m_P2;};

	///Get the distance ot P1() and P2().
	double distance()const{return m_dist;};

	///Get the tangent as P1().
	const MGUnit_vector& tan1()const{return m_tan1;};

	///Get the tangent as P2().
	const MGUnit_vector& tan2()const{return m_tan2;};

	///tan1 and tan2's angle in radian.
	double tandiff()const{return m_tandiff;};

	///Get curvature direction at P1().
	const MGUnit_vector& normal1()const{return m_normal1;};

	///Get curvature direction at P2().
	const MGUnit_vector& normal2()const{return m_normal2;};
	
	///Get the normal1 and normal2's angle in radian.
	double normaldiff()const{return m_normaldiff;};

	///Get curvature at P1().
	double curvature1()const{return m_curvature1;};

	///Get curvature at P2().
	double curvature2()const{return m_curvature2;};

private:
	/// continuity indicator
	CONTINUITY m_continuity; ///< geometric continuity

	double m_param1; ///<parameter value of start or end of curve1.
	double m_param2; ///<parameter value of start or end of curve2.

	MGPosition m_P1; ///<positional data at m_param1.
	MGPosition m_P2; ///<positional data at m_param2.
	double m_dist; ///< distance between two curve's endpoints.

	MGUnit_vector m_tan1; ///< tangent at m_param1
	MGUnit_vector m_tan2; ///< tangent at m_param2
	double m_tandiff; ///< angle between [t1] and [t2] in radian

	double m_curvature1; ///< curvature at m_param1
	double m_curvature2; ///< curvature at m_param2
	MGUnit_vector m_normal1; ///< curvature direction at m_param1
	MGUnit_vector m_normal2; ///< curvature direction at m_param2
	double m_normaldiff;///< angle between the two curvature direction in radian
};

/** @} */ // end of GEORelated group
#endif //__MGCurveContinuity_H__
