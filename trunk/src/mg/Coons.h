/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCoons_HH_
#define _MGCoons_HH_

#include <iosfwd>
#include "mg/Pvector.h"

// MGCoons.h

// Forward Declaration
class MGCurve;
class MGLBRep;
class MGSPointSeq;

/** @addtogroup GEORelated
 *  @{
 */

/// Defines Coons Patch surface.
class MGCLASS MGCoons {
 
public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGCoons& );

//////////// Constructor ////////////

//Default constructor.
//MGCoons(){;};

MGCoons(
	MGPvector<MGLBRep>& perimeters,
	MGPvector<MGLBRep>& derivatives
);

MGCoons(
	MGPvector<MGCurve>& perimeters,
	MGPvector<MGCurve>& derivatives
);

//Copy constructor.
//MGCoons(const MGCoons& rhs);

////////////Destructor/////////

//~MGCoons();
									
//////////// Operator Overload ////////////

//MGCoons& operator=(const MGCoons&); ///Assignment

//////////// Member Function ////////////

///get the derivative data d2f/((dv)(du)) at corner i.
const MGVector& d2fdvu(size_t i)const;

///Return perimeter i's derivative data.
const MGCurve& derivative(size_t i)const{return *(m_derivatives[i]);};

///Evaluate surface data.
///Currently ndu=ndv=0 is assumed.
MGVector eval(
	double u, double v	///< Parameter value of the surface,
						///< must be 0<=u,v<=1.
	, size_t ndu=0		///< Order of derivative along u.
	, size_t ndv=0		///< Order of derivative along v.
	) const;

void eval(
	const MGNDDArray&	utau,		///<u方向のデータポイント
	const MGNDDArray&	vtau,		///<v方向のデータポイント
	MGSPointSeq&		spoint///<evaluated data will be output to spoint.
)const;

///Get space dimension.
size_t sdim()const;

//////////// Member Data ////////////

private:

	MGPvector<MGCurve> m_perimeters;///<perimeters.
		///<m_perimeters[0]:v-min, [1]:u-max, [2]:v-max, [3]:u-min.
		///<the parameter ranges of all the curves must be [0,1].
	MGPvector<MGCurve> m_derivatives;///<Derivatives along perimeters.
		///<m_derivatives[0]:along perimeters[0](df/dv(u,0))
		///<[1]:along perimeters[1](df/du(1,v))
		///<[2]:along perimeters[2](df/dv(u,1))
		///<[3]:along perimeters[3](df/du(0,v))

	///Data at each corner
	MGVector m_f00, m_f01, m_f10, m_f11;
	MGVector m_dfdu00, m_dfdu01, m_dfdu10, m_dfdu11;
	MGVector m_dfdv00, m_dfdv01, m_dfdv10, m_dfdv11;
	MGVector m_df2duv00, m_df2dvu00;
	MGVector m_df2duv01, m_df2dvu01;
	MGVector m_df2duv10, m_df2dvu10;
	MGVector m_df2duv11, m_df2dvu11;

void eval_corner();
MGVector eval_inner(
	double u, double v	///< Parameter value of the surface.
						///< must be 0<=u,v<=1.
	) const;

friend double mgHermite0(double t);
friend double mgHermite1(double t);
friend double mgHermite2(double t);
friend double mgHermite3(double t);

};

/** @} */ // end of GEORelated group
#endif
