/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSBRepVecTP_HH_
#define _MGSBRepVecTP_HH_

#include <assert.h>
#include "mg/Pvector.h"
#include "mg/LBRep.h"
#include "mg/SBRepTP.h"

// MGSBRepVecTP.h
//

class MGSurface;
class MGLBRep;
class MGOfstream;
class MGIfstream;

/** @addtogroup GEORelated
 *  @{
 */

///Defines Tangent Plane Line B-Representation Class.
///Tangent plane is a line b-representation of (unit)normal vector of
///tangent plane along a surface perimeter.
///MGSBRepVecTP has a vector of 4 newed object pointer ofMGLBRep representing the normal vetor
///B-rep along the 4 perimeter of surface. And regarding to their ownership
///MGSBRepVecTP acts just like std::auto_ptr. That is, all of the ownership will be
///transfered to the copied, or assigned object.
///See the copy oroperator=() below.
class MGCLASS MGSBRepVecTP{

public:

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream& ostrm, const MGSBRepVecTP& vectp);

//////////// Constructor ////////////

///Default Constructor, will be set as no TPs' are specified.
MGSBRepVecTP(){;};

///Copy Constructor.
///all of the ownership of m_TP of tp2 will be
///transfered to this object.
MGSBRepVecTP(const MGSBRepVecTP& vectp2);

///Conversion constructor from ordinary MGSBRepTP.
MGSBRepVecTP(const MGSBRepTP& tp2);
	  
//////////// Destructor ////////////

//~MGSBRepVecTP();

//////////// Operator overload. ////////////

///Assignment.
///all of the ownership of m_TP of tp2 will be
///transfered to this object.
MGSBRepVecTP& operator=(const MGSBRepVecTP& vectp2);

//////////// Member Function ////////////

///change the parameter range to (t0,t1).
void change_range(
	bool along_u,	///<the objective range is along u parameter or v.
	double t0, double t1
);
void change_range(
	bool along_u,	///<the objective range is along u parameter or v.
	const MGInterval& prange
);

///evaluate TP at the perimeter i's parameter t.
///Function's return value is:
///true if t was inside the parameter range of a tangent plane of m_TP[i].
///false if t was outside the parameter range of the m_TP[i] for all i.
bool eval(
	size_t i,		///<perimeter numeber.
	double t,		///<parameter vaule of perimeter i.
	MGVector& normal///<evaluated normal will be returned.
)const;

/// Compute the maximum (absolute) cos value of between vector deris[i](t) 
/// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
/// parameter of the data point obtained from deris[i]'s knot vector.
///Function's return value is the max out of cosmax[.].
double get_perimeters_max_cos(
	const MGPvector<MGLBRep>& deris,     ///< the size must be 4
	double                    taumax[4], ///< parameter on which the maximum value attains will be stored
	double                    cosmax[4]  ///< the maximum value will be stored
)const;

/// Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
/// and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
/// the same as srf.perimeter_curve(i), and t is a common parameter
/// of deris[i] and TP(i).
///Function's return value is the max out of sinmax[.].
double get_perimeters_max_sin(
	const MGSurface& srf,       ///< surface which must corresponds to this object
	double         taumax[4], ///< parameters on which the maximum value attains will be stored.
	double         sinmax[4],  ///< the maximum value will be stored.
	bool*          eval=0	///<indicates perimeters to evalate if eval!=null
			///<When eval[i] is true, perimeter i is evaluated for 0<=i<=3.
)const;

///Return true if at least one TP is specified at i-th perimeter.
///i=0, 2 are v=min and max u-parameter line.
///i=1, 3 are u=max and min v-parameter line.
bool specified(size_t i) const{assert(i<4); return m_TP[i].size()>0;};

///Set i-th perimeter's TP as a null, as an unspecified one.
void set_TP_null(size_t i);

///Set i-th perimeter's TP(copy version).
///vectp[i] must be newed objects, and all of the ownership will be transferer to
///this instance.
void set_TP(
	size_t i,					///<perimeter numeber.
	MGPvector<MGLBRep>& vectp,
	const MGInterval& prange		///<Whole perimeter's parameter range.
);

///Set i-th perimeter's TP(std::vector version).
///vectp[i] must be newed objects, and all of the ownership will be transferer to
///this instance.
void set_TP(
	size_t i,					///<perimeter number.
	std::vector<MGLBRep*>& vectp,
	const MGInterval& prange		///<Whole perimeter's parameter range.
);

///Return i-th perimeter's TP.
const MGPvector<MGLBRep>& vecTP(size_t i) const{assert(i<4);return m_TP[i];}
MGPvector<MGLBRep>& vecTP(size_t i) {assert(i<4);return m_TP[i];}

MGPvector<MGLBRep>* vecTP(){return m_TP;};

private:
//////////// Member Data ////////////
	MGPvector<MGLBRep> m_TP[4];	///<Tangent Plane will be stored,
				///<Tangent plane m_TP[i] is a vector of line b-representations of
				///<(unit)normal vector along the i-th perimeter,
				///<Parameter range of the TP is the same as u or v parameter
				///<range of the corresponding surface representation,
	///< m_TP[0]: v=min boundary line, m_TP[1]: u=max boundary line,
	///< m_TP[2]: v=max boundary line, m_TP[3]: u=min boundary line.

	MGInterval m_prange[4];///m_prange[j] is the parameter range of perimeter j

	bool m_to_SE[8];///m_to_SE[i*2] and [i*2+1] indicate if m_TP[i]'s parameter range covers
		///from the start point(i*2), or to the end point(i*2+1),
		///when m_to_SE[i*2] is true, m_TP[i] starts from the start point of perimeter i,
		///when m_to_SE[i*2+1] is true, m_TP[i] ends at the end point of the perimeter i,
		///when false, m_TP[i] starts from the inner point, or ends at the inner point. 
};

/** @} */ // end of GEORelated group
#endif
