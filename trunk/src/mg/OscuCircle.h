/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGOscuCircle_HH_
#define _MGOscuCircle_HH_
/** @addtogroup GEORelated
 *  @{
 */

#include "mg/OscuCircleData.h"
#include <vector>

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGOscuCircleData>;
#pragma warning( pop )
#endif

// MGOscuCircle.h
//

//Forward Declaration
class MGIfstream;
class MGOfstream;

/// Defines Array of OscuCircle data.
///This class is used for MGLBRep constructor:
///	MGLBRep(						///BLGCS
///			const MGLBRepEndC& begin,	///Begin end condition
///			const MGLBRepEndC& end,		///End end conditoion
///			const MGBPointSeq& points,	///Point seq data
///			const int* point_kind,		///Point kind of above point.
///			const MGOscuCircle& circle,	///Provides osculating circle data.
///			int &error);				///Error flag.
///Defines array of MGOscuCircleData which provides index of points and
///radius of the osculating circles.
class MGCLASS MGOscuCircle {

public:

//////////// Constructor ////////////

/// Dummy constructor, setts m_n=0
MGOscuCircle():m_n(0){;};

///	MGOscuCircle(const MGOscuCircle&);	///Copy Constructor.

//////////// Destructor ////////////
//	~MGOscuCircle(){;};	

//////////// Operator overload. ////////////

//	MGOscuCircle& operator =(MGOscuCircle&);///Assignment operator difinition.
//                        We use default opassignment operator overload 

///Reference i-th osculating circle data.
const MGOscuCircleData& operator ()(size_t i) const;

//////////// Member Function ////////////

///Add to the end of list.
MGOscuCircle& add(const MGOscuCircleData&);	

///Add to the end of list.
MGOscuCircle& add(size_t index, double radious);

///Remove i-th OscuCircleData.
MGOscuCircleData remove(size_t i);

size_t length() const {return m_n;};

///Debug Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGOscuCircle& );

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

//////////// Member Data ////////////

private:
	size_t m_n;			///<Number of stored OscuCircle data.
	std::vector<MGOscuCircleData> m_circle;	///<OscuCircle data list.

};

/** @} */ // end of GEORelated group
#endif
