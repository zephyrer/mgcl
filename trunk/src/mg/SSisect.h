/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSSisect_HH_
#define _MGSSisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Curve.h"

// MGSSisect.h
// Header for MGSSisect

//Forward Declaration
class MGCurve;

///MGSSisect represents one intersection line of two surfaces.
///A list of intersection lines are expressed by MGSSisect_list.
///The behavior of MGSSisect is like a auto_ptr. Copy or assignment
///of MGSSisect means transfer of the ownership of all the included curve
///to copied or assigned MGSSisect and original MGSSisect does not have the
///ownership of the curves any more. User should be aware of it.
/// Surface と Surface の交線を一つのみ表現する。交線の集合は別に表現される。
class MGCLASS MGSSisect:public MGisect{

public:

//////////// Constructor ////////////

///Void constructou. 初期化なしでDummy交線を生成
MGSSisect()
:m_iline(0), m_param1(0), m_param2(0),m_rel(MGSSREL_UNKNOWN)
{;};

/// Copy Constructor;
/// ssi's ownership of all the three curves will be released.
MGSSisect(const MGSSisect& ssi);

///Construct providing all the raw data.
///The ownership of iline, param1, and param2 are all transfered to MGSSisect.
///All of these objects must be newed ones.
MGSSisect(
	MGCurve* iline,	///<Pointer of newed object.
	MGCurve* param1,///<Pointer of newed object.
	MGCurve* param2,///<Pointer of newed object.
	const MGSSRELATION r1=MGSSREL_UNKNOWN)
	:m_iline(iline), m_param1(param1), m_param2(param2), m_rel(r1){;};

///Construct providing all the raw data.
///Copy version. Copy of the three curves will take place.
MGSSisect(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1=MGSSREL_UNKNOWN);

//////////// Destructor ////////////
~MGSSisect();

//////////// Operator overload ////////////

///Assignment
/// ssi's ownership of all the three curves will be released.
MGSSisect& operator= (const MGSSisect& ssi);

///Comparison operator.
bool operator< (const MGSSisect& ssi2)const;
bool operator> (const MGSSisect& ssi2)const{return ssi2<(*this);};
bool operator<= (const MGSSisect& ssi2)const{return !(ssi2<(*this));};
bool operator>= (const MGSSisect& ssi2)const{return !((*this)<ssi2);};
bool operator== (const MGSSisect& ssi2)const;
bool operator!= (const MGSSisect& ssi2)const{return !operator==(ssi2);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
bool operator< (const MGCFisect& is)const{return false;};
bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

//////////// Memeber Function ////////////

///Test if two ssi's world curve have common parts (in line_zero()).
///Fucntion's return value is 
///		1:have common part.
///		0>=:no common part(except a point).
int has_common(const MGSSisect& ssi2)const;

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12(){replace12();};

///Test if this SSI is null.
bool is_null()const{return m_iline==0;};

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return *m_iline;};

///Return the 1st object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
const MGCurve* isect1_param1()const{return m_param1;};

///Return the 2nd object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
const MGCurve* isect1_param2()const{return m_param2;};

/// 交線の座標値表現を返却する
///Return (x,y,z) coordinate representation intersection line.
MGCurve& line() const{return *m_iline;}

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 1;};

///negate the direction of the intersection line.
void negate();

/// Output function.
std::ostream& out(std::ostream& ostrm)const;

/// 交線の第１surface のパラメータ値表現を返却する
///Return (u,v) parameter representation intersection line
///of the 1st surface.
MGCurve& param1() const{return *m_param1;}

/// 交線の第２surface のパラメータ値表現を返却する
///Return (u,v) parameter representation intersection line
///of the 2nd surface.
MGCurve& param2() const{return *m_param2;}

/// 交線での両surface の関係を返却。
///Return the relationship at the intersection line.
MGSSRELATION rel() const{return m_rel;}

///Release each curve pointer from this.
///After the use of release_xxxx(), MGSSisect does not have the ownership of
///the each curve.
MGCurve* release_line();
MGCurve* release_param1();
MGCurve* release_param2();

///Replace 1st and 2nd order of the parameter line representation.
MGSSisect& replace12();

void set_null();

private:

//////////// Member Data ////////////
	
	MGCurve* m_iline;	///<(x,y,z)coordinate representaion of the line.
					///< New-ed curve will be stored.
					///<交線の座標値よる表現
	MGCurve* m_param1;	///<(u,v) representaion of the line of the first
					///< New-ed curve will be stored.
					///<surface. 交線の第１surface のパラメータ(u,v)による表現
	MGCurve* m_param2;	///<(u,v) representaion of the line of the second
					///< New-ed curve will be stored.
					///<surface. 交線の第2surface のパラメータ(u,v)による表現
	MGSSRELATION m_rel;	///<Two surfaces relationship at the intersection line.
					///<交線における両Surfaceの関係。

	friend class MGSSisect_list;
};

/** @} */ // end of IsectContainer group
#endif
