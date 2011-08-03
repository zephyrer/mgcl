/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGHHisect_HH_
#define _MGHHisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/CompositeCurve.h"
#include "mg/FPline.h"

#if defined(MGCL_DLL)
#include "mg/DequeProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGDequeProxy<MGFPline>;
#pragma warning( pop )
#else
#include <deque>
#endif

class MGSSisect;
class MGFSurface;
class MGHHisect_vector;

///MGHHisect is to represent one continuous intersection line of shells.
///Intersection lines a shell with a shell, a face, or a surface.
///(MGCompositeCurve* iline, deque<MGFPline> uvl1, deque<MGFPline> uvl2)
///where iline is a world coordinate rep of the line, uvl1 is a deque of
///1st shell's face parameter rep, uvl2 is a deque of the 2nd shell's or
///face's parameter rep of the intersection line.
///uvl1[i] corresponds to uvl2[i] one by one for all i.
///The parameter ranges of all the uvl1[i] are continuous and the total of them is
///equal to the parameter range of iline. For uvl2, the same.
///Let uvl1[i]'s start parameter be t1, and end parameter t2, then
///uvl2[i]'s parameter range is also from t1 to t2.
///Let sf1 be MGSurfCurve(f1's surface, uvl1[i]), then sf1 is the same curve as
///the iline's part of the parameter range t1 to t2. And sf1 is also equal to
///MGSurfCurve(f2's surface, uvl2[i]).
///MGHHisect uses MGFPline to represent the intersection lines.
///The behavior of MGHHisect(and MGFPline) is like an auto_ptr. Copy or assignment
///of MGHHisect means transfer of the ownership of all the included curves
///to copied or assigned MGHHisect and original MGHHisect does not have the
///ownership of the curves any more. Users should be aware of it.
///MGHHisect is also used to represent a pojection curve. In this case the size of
///uvline2 is zero.
///**** Projection line rep and intersection line rep cannot be mixed. ****
class MGCLASS MGHHisect{

public:
#if defined(MGCL_DLL)
	typedef MGDequeProxy<MGFPline> container_type;
#else
	typedef std::deque<MGFPline> container_type;
#endif

	typedef container_type::iterator iterator;
	typedef container_type::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream& ostrm, const MGHHisect& hhi);
	
/////////Constructor/////////
MGHHisect():m_iline(0){;};	///void constructor.

///Copy constructor.
MGHHisect(const MGHHisect& hhi);

///Construct from MGSSisect.
MGHHisect(
	const MGFSurface* face1,	///<face1. This must not be null.
	const MGFSurface* face2,	///<face2. This may be null
							///<(e.g. for face2 that is actually a surface).
	MGSSisect& ssi);		///<intersection line of face1 and face2 expressed as
							///<MGSSisect.

///Construct two faces intersection lines.
///uvline1 and 2 makes uvlines whose vector length is 1.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
///When face2 is null, and uvlines2!=null, it indicates face2 is actually a surface.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
MGHHisect(
	MGCurve* iline,		///<Intersection line of world coordinates.
			///<iline must be a newed object pointer and MGHHisect takes the ownership.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,	///<Intersection line of face1's (u,v) coordinates.
			///<uvline1 must be a newed object pointer and MGHHisect takes the ownership.
	const MGFSurface* face2=0,///<When face2 is null, and uvlines2!=null, it indicates
						///<face2 is actually a surface.
	MGCurve* uvline2=0  ///<Intersection line of face2's (u,v) coordinates.
						///<Takes the ownership.
);

/////////Destructor/////////
~MGHHisect();

/////////Operator overload/////////

///Assignment.
MGHHisect& operator=(const MGHHisect& hhi);

///Comparison operator.
bool operator< (const MGHHisect& hhi2)const;
bool operator> (const MGHHisect& hhi2)const{return hhi2<(*this);};
bool operator<= (const MGHHisect& hhi2)const{return !(hhi2<(*this));};
bool operator>= (const MGHHisect& hhi2)const{return !((*this)<hhi2);};
bool operator== (const MGHHisect& hhi2)const;
bool operator!= (const MGHHisect& hhi2)const{return !operator==(hhi2);};

/////////Member function/////////

///Extract a connected line from hhivec one by one and build one continuous
///line. Extracted lines will be released from hhivec.
void build_one(MGHHisect& hhi);
void build_one(MGHHisect_vector& hhivec);

///Change parameter range, be able to change the direction by providing
///t0 greater than t1.
void change_range(
	double t0,		///<Parameter value for the start of original. 
	double t1);	///<Parameter value for the end of original. 

///Connect a line to this HHisect.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
///****Projection lines rep. and intersection lines rep cannot be mixed.
///iline, uvline1, and uvline2(or hhi2) must have the same direction.
///iline's direction must be equal to this HHisect's.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
void connect_line_to_end(
	MGCurve* iline,		///<Intersection line of world coordinates.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,	///<Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2=0,///<face2. This may be null
						///<(e.g. for face2 that is actually a surface).
	MGCurve* uvline2=0///<Intersection line of face2's (u,v) coordinates.
);						///<takes the ownership of all the curves of ssi.
void connect_line_to_end(
	MGHHisect& hhi2		///<After connected, this hhi2's member data's ownership
						///<will be transfered to this MGHHisect,
						///<just like std::auto_ptr's assignment.
);

///Connect a line to this HHisect.
///When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
///used. This case occurs when MGHHisect is used to represent projection lines.
///****Projection lines rep. and intersection lines rep cannot be mixed.
///iline, uvline1, and uvline2(or hhi2) must have the same direction.
///iline's direction must be opposite to this HHisect's.
///MGHHisect takes the ownership of iline, uvline1, and uvline2
///(These must be newed objects).
void connect_line_to_start(
	MGCurve* iline,		///<Intersection line of world coordinates.
	const MGFSurface* face1,///<face1. This must not be null.
	MGCurve* uvline1,	///<Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2=0,///<face2. This may be null
						///<(e.g. for face2 that is actually a surface).
	MGCurve* uvline2=0);///<Intersection line of face2's (u,v) coordinates.
						///<takes the ownership of all the curves of ssi.
void connect_line_to_start(
	MGHHisect& hhi2		///<After connected, this hhi2's member data's ownership
						///<will be transfered to this MGHHisect,
						///<just like std::auto_ptr's assignment.
);

///Test if this has face2's information.
bool has_face2_data()const{return !m_uvlines2.empty();};

///Return the world coordinate isect data.
const MGCurve& iline()const{return *m_iline;};
MGCurve& iline(){return *m_iline;};

///Test if this is a null HHisect or not.
bool is_null()const{return m_iline==0;};

double param_e()const{return m_iline->param_e();};
double param_s()const{return m_iline->param_s();};

///Reverse the direction of this intersection line.
void reverse_direction();

///return number of uvlines.
size_t num_of_uvline() const{return m_uvlines1.size();};

///Release the pointer of the last curve.
///Returned will be the released MGCurve pointer.
void release_back(
	MGCurve*& ilineLast,
	MGFPline uvline1Last,
	MGFPline uvline2Last
);

///Release the pointer of the 1st curve.
///Returned will be the released MGCurve pointer.
void release_front(
	MGCurve*& iline1st,
	MGFPline uvline11st,
	MGFPline uvline21st
);

///Release the pointer of the iline curve.
///Returned will be the released MGCurve pointer.
MGCompositeCurve* release_line(){MGCompositeCurve* a=m_iline; m_iline=0; return a;};

///Return uvline1.
const std::deque<MGFPline>& uvlines1() const{
#if defined(MGCL_DLL)
	return m_uvlines1.std_container();
#else
	return m_uvlines1;
#endif
}

///Return uvline2.
const std::deque<MGFPline>& uvlines2() const{
#if defined(MGCL_DLL)
	return m_uvlines2.std_container();
#else
	return m_uvlines2;
#endif
}

///Return i-th uvline.
const MGFPline& uvline1(size_t i)const;
const MGFPline& uvline2(size_t i)const;

///exchange12 1st and 2nd lines.
///This can be used only for intersection line rep.
void exchange12();

private:
	mutable MGCompositeCurve* m_iline;
			///<World coordinates (x,y,z) rep of the intersection line.
	container_type m_uvlines1;
			///<1st (u,v) parameter rep line of the intersection line.
	container_type m_uvlines2;
			///<2nd (u,v) parameter rep line of the intersection line.
			///<This vector lenght can be 0 even m_uvlines1's vector length is more than 0,
			///<which means this isect does not have 2nd face data.

};

/** @} */ // end of IsectContainer group
#endif
