/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFPline_HH_
#define _MGFPline_HH_

class MGCurve;
class MGFSurface;

/** @addtogroup GEORelated
 *  @{
 */

///Define MGFPline Class, Face's (u,v) parameter value line.
///MGFPline is to represent an parameter (u,v) line of a face.
///(MGFSurface* f, MGCurve uvline) where f is a face pointer, and uvline is
///a parameter (u,v) line of the face f.
///
///MGFPline is used to express Shell's intersection lines.
///The behavior of MGFPline is like an auto_ptr. Copy or assignment
///of MGFPline means transfer of the ownership of all the included curves
///to copied or assigned MGFPline and original MGFPline does not have the
///ownership of the curves any more. Users should be aware of it.
class MGCLASS MGFPline{

public:
	
///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream& ostrm, const MGFPline& fpl);

////////Constructor/////////
MGFPline():m_face(0),m_uvline(0){;};	///void constructor.

///Construct from all the necessary data.
MGFPline(
	const MGFSurface* face,	///face1.
	MGCurve* uvline)///(u,v) line of the face, takes the ownership of the curve.
					///That is uvline must be a newed object pointer.
	:m_face(face), m_uvline(uvline){;};

/// Copy Constructor;
/// fpl's ownership of all the curve will be transfered to
/// the new MGFPline object.
MGFPline(const MGFPline& fpl);

//////////// Destructor ////////////
~MGFPline();	///uvline will be deleted.

////////Operator oveload////////

///Assingment.
///The ownership of the curve in fpl will be transfered to this MGFPline.
MGFPline& operator= (const MGFPline& fpl);

///Comparison operator.
bool operator< (const MGFPline& fpl2)const;
bool operator> (const MGFPline& fpl2)const{return fpl2<(*this);};
bool operator<= (const MGFPline& fpl2)const{return !(fpl2<(*this));};
bool operator>= (const MGFPline& fpl2)const{return !((*this)<fpl2);};
bool operator== (const MGFPline& fpl2)const;
bool operator!= (const MGFPline& fpl2)const{return !operator==(fpl2);};

////////Member function////////

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t0,	///<Parameter value for the start of original. 
	double t1	///<Parameter value for the end of original. 
);

///Get face's pointer.
const MGFSurface* face()const{return m_face;};

///Reverse the direction of this line.
void reverse_direction();

///Release the uvline curve pointer from this.
///After the use of release_line(), MGFPline does not have the ownership of
///the curve.
MGCurve* release_line();

///Return face's (u,v) parameter representation line.
const MGCurve& uvline() const{return *m_uvline;}
MGCurve& uvline() {return *m_uvline;}

private:
	const MGFSurface* m_face;	///<Face pointer.
	mutable MGCurve* m_uvline;///<2D line whose coordinates are (u,v) of the m_face.
						///<m_uvlline is a newed object, and the ownership is 
						///<controled by this class.

};

/** @} */ // end of GEORelated group

#endif
