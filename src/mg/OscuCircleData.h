/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGOscuCircleData_HH_
#define _MGOscuCircleData_HH_
/** @addtogroup GEORelated
 *  @{
 */
#include "mg/MGCL.h"

// MGOscuCircleData.h
//

class MGIfstream;
class MGOfstream;

///The class for MGLBRep constructor of osculating circles.
///This class is used for MGLBRep constructor:
///	MGLBRep(						///BLGCS
///			const MGLBRepEndC& begin,	///Begin end condition
///			const MGLBRepEndC& end,		///End end conditoion
///			const MGBPointSeq& points,	///Point seq data
///			const int* point_kind,		///Point kind of above point.
///			const MGOscuCircle& circle,	///Provides osculating circle data.
///			int &error);				///Error flag.
/// Defines OscuCircle data, index of BPointSeq points and circle radius.
class MGCLASS MGOscuCircleData {

public:

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGOscuCircleData& );

////////////Constructor////////////

///Default Constructor
MGOscuCircleData(){;};

MGOscuCircleData( size_t index,		///Index of BPointSeq indicating where
									///circle be inserted.
				double radius);		///Radius of circle.
//Destructor
//	~MGOscuCircleData();	

///Operator overload.

bool operator< (const MGOscuCircleData& ocd2)const;
bool operator> (const MGOscuCircleData& ocd2)const{return ocd2<(*this);};
bool operator<= (const MGOscuCircleData& ocd2)const{return !(ocd2<(*this));};
bool operator>= (const MGOscuCircleData& ocd2)const{return !((*this)<ocd2);};
bool operator== (const MGOscuCircleData& ocd2)const;
bool operator!= (const MGOscuCircleData& ocd2)const{return !operator==(ocd2);};

////////////Member Function////////////////
	
size_t index() const {return m_index;}
double radius() const {return m_radius;}

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

////////////Member Data//////////
private:
	size_t m_index;			///<Index of points.
	double	m_radius;		///<radius of oscurating circle.

};

/** @} */ // end of GEORelated group
#endif
