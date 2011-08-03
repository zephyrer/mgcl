/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFPoint_HH_
#define _MGFPoint_HH_

#include "mg/Position.h"
class MGFSurface;

/** @addtogroup GEORelated
 *  @{
 */

///MGFPoint is to represent a Face or Surface point.
///The expression is:
///(MGFace* f, MGPosition uv), where f is a face pointer of interest, 
///and uv is the parameter value of the face f.
class MGCLASS MGFPoint{

public:
	
///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGFPoint& );

////////Constructor/////////
MGFPoint():m_face(0){;};	///void constructor.

///Construct from all the necessary data.
MGFPoint(
	const MGFSurface& face,	///<face.
	const MGPosition& uv	///<Parameter values of the face.
	):m_face(&face), m_uv(uv){;};

////////Operator oveload////////

///The sequence of MGFPoint's are defined in the same face's MGFPoint.
///If two MGFPoint do not belong to the same face, their sequences are undefined.
///When FP1(u1,v1) and FP2(u2, v2) are two MGFPoint, then FP1<FP2
///if u1<u2, or (u1=u2 and v1<v2).
bool operator< (const MGFPoint& fp)const;
bool operator> (const MGFPoint& fp)const{return fp<(*this);};
bool operator<= (const MGFPoint& fp)const{return !(fp<(*this));};
bool operator>= (const MGFPoint& fp)const{return !((*this)<fp);};

///Two Fpoints are equal if they belong to one face and their distance in parameter
///(u1,v2) and (u2,v2) is less than parameter_error() of the face.
bool operator== (const MGFPoint& fp)const;
bool operator!= (const MGFPoint& fp)const{return !operator==(fp);};

////////Member function////////

///Evaluation of the Face at the FPoint.
///When nderi=0, get the positional data at the point.
MGVector eval(size_t ndu=0, size_t ndv=0)const;

///return the face.
const MGFSurface& fsurface()const{return *m_face;};

///Return the parameter value of the face.
const MGPosition& uv()const{return m_uv;};

private:
	const MGFSurface* m_face;///MGFSurface pointer of the point (u,v).
	MGPosition m_uv;		///m_face's parameter value (u,v).

};

//////////// GLOBAL FUNCTIONS FOR SHELL /////////////////////

///Evaluate face's (shell's) point from the MGFPoint.
MGDECL MGVector eval(
	const MGFPoint& uv,	///<Face point.
	size_t ndu=0, size_t ndv=0	///<Order of derivative along u and v direction.
);

///Compute normal vector(not unit) at uv.
MGDECL MGVector normal(const MGFPoint& uv);

///Compute unit normal vector at uv.
MGDECL MGUnit_vector unit_normal(const MGFPoint& uv);

/** @} */ // end of GEORelated group

#endif
