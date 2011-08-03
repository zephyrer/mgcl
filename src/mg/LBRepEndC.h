/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLBRepEndC_HH_
#define _MGLBRepEndC_HH_
#include "mg/MGCL.h"
#include "mg/Vector.h"

// MGLBRepEndC.h
//

class MGNDDArray;
class MGBPointSeq;
class MGCurve;
class MGIfstream;
class MGOfstream;

/** @addtogroup GEORelated
 *  @{
 */

/// Defines End Condition of Line B-Representation.
///End condition to get spline by interpolation.
///enum MGENDCOND {
///	MGENDC_UNKNOWN=0,	/// Unknown(usually not used).–¢’m
///	MGENDC_1D  =1,		/// 1st deravative provided.
///	MGENDC_2D  =2,		/// 2nd deravative provided.
///	MGENDC_NO  =3,		/// no end cond(only positional data)
///	MGENDC_12D =4		/// both 1st and 2nd deravatives provided.
///};
class MGCLASS MGLBRepEndC {

public:

////////////Constructor////////////

///Default Constructor.
MGLBRepEndC():m_cond(MGENDC_NO){;};

///Construct of m_cond=MGENDC_1D or MGENDC_2D.
MGLBRepEndC(
	MGENDCOND cond,	///<Type of end condition
					///<(MGENDC_1D or MGENDC_2D).
	const MGVector& deriv	///<Derivative inf according to cond
);

///Construct of m_cond=MGENDC_12D
MGLBRepEndC(
	const MGVector& first_deriv,	///<1st derivative
	const MGVector& second_deriv	///<2nd derivative
);

/// Given data points ordinates and abscissa, compute approximate
/// 1st derivative.
MGLBRepEndC(
	int start,			///<Indicates start(start==true)
						///< condition or end.
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	int &error				///<Error flag.
);

/// Given Positional data sequence, compute approximate
/// 1st derivative. The 1st derivative is unit vector.
MGLBRepEndC(
	int start,			///<Indicates start(start==true)
						///< condition or end.
	const MGBPointSeq& points	///<Point seq data
);

/// Given MGCurve, construct the curve's end condition.
MGLBRepEndC(
	int start,		///<Indicates start(start==true) condition or end.
	MGENDCOND cond,	///<Type of end condition(MGENDC_1D, MGENDC_2D, or MGENDC_12D)
	const MGCurve& curve///<Curve
);

///Destructor
///	~MGLBRepEndC();	  We use default constructor.

///Debug Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGLBRepEndC& );

////////////Member Function////////////

MGENDCOND cond() const {return m_cond;}
const MGVector& first() const {return m_1deriv;}
MGVector& first() {return m_1deriv;}
const MGVector& second() const {return m_2deriv;}
MGVector& second() {return m_2deriv;}

///Initialize the instance. Will be set to the same as constructed by the void 
///constructor.
void initialize();

///Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
void set_1st(const MGVector& first_deriv);

///Set 2nd deriv and change condition type to MGENDC_2D or MGENDC_12D.
void set_2nd(const MGVector& second_deriv);

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Member Data
private:

	MGENDCOND m_cond;	///<Type of end condition.
	MGVector  m_1deriv;	///<1st derivative stored 
						///<when m_cond=MGENDC_1D or MGENDC_12D
	MGVector  m_2deriv;	///<2nd derivative stored
						///<when m_cond=MGENDC_2D	or MGENDC_12D

};

/** @} */ // end of GEORelated group

#endif
