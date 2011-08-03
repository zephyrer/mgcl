/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCCisect_HH_
#define _MGCCisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Point.h"

// MGCCisect.h
// Header for MGCCisect

///CCisect is for one intersection representation of two curves.
///Curve �� Curve �̌�_����_�̂ݕ\������.
///��_�̏W���͕ʂɕ\�������B
///�{�N���X�͒��� �� ���� �̌�_�ȂǁA���X��_�̌�_�̕ԋp�p�ɗ��p�����B
///If more than one are necessary to hold, CCisect_list should be used.
class MGCLASS MGCCisect:public MGisect{

public:

////////////Constructor////////////

///Void Constructor(�������Ȃ��Ō�_�𐶐�)
MGCCisect();

///Input all necessary components(�S�ẴR���|�[�l���g���w�肵�Č�_�𐶐�)
MGCCisect(const MGPosition &is, double t1, double t2,
		  const MGCCRELATION r1=MGCCREL_UNKNOWN);

////////////Operator Overload////////////

bool operator< (const MGCCisect& cci)const{return m_param1<cci.m_param1;};
bool operator> (const MGCCisect& cci)const{return m_param1>cci.m_param1;};
bool operator<= (const MGCCisect& cci)const{return m_param1<=cci.m_param1;};
bool operator>= (const MGCCisect& cci)const{return m_param1>=cci.m_param1;};
bool operator== (const MGCCisect& cci)const;
bool operator!= (const MGCCisect& cci)const{return !operator==(cci);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCSisect& is)const{return true;};
bool operator< (const MGCFisect& is)const{return true;};
bool operator< (const MGSSisect& is)const{return true;};
bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

///Debug Function
/// Output virtual function.
std::ostream& out(std::ostream& ostrm)const;

//////////////Memeber Function////////////

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12();

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return m_ipoint;};

///Return the 1st object's parameter value of the intersection.
MGPosition isect0_param1()const{return MGPosition(1,&m_param1);};

///Return the 2nd object's parameter value of the intersection.
MGPosition isect0_param2()const{return MGPosition(1,&m_param2);};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 0;};

///Return Two curves' relationship(��_�ł̗�curve �̊֌W��ԋp)
MGCCRELATION rel() const{return m_rel;}

///Return coordinate values(��_�̍��W�l��ԋp����)
const MGPosition& point() const{return m_ipoint.position();}

///Return parameter value of 1st curve.
///��_�̑�Pcurve �̃p�����[�^�l��ԋp����
double param1() const{return m_param1;}

///Return parameter value of 2nd curve.
/// ��_�̑�Qcurve �̃p�����[�^�l��ԋp����
double param2() const{return m_param2;}

///Set param1 data.
void set_param1(double t1){ m_param1=t1;}

///Set param2 data.
void set_param2(double t2){ m_param1=t2;}

private:
///Member Data
	MGPoint m_ipoint; ///< coordinate values(��_�̍��W�l)
	double m_param1;	 ///< parameter value of 1st curve
						 ///<(��_�ɂ������Pcurve �̃p�����[�^�l)
	double m_param2;	 ///< parameter value of 2nd curve
						 ///< ��_�ɂ������Qcurve �̃p�����[�^�l
	MGCCRELATION m_rel;	 ///< Two curves' relationship at the i.p.,
						 ///<��_�ɂ����闼�J�[�u�̊֌W�B

};

/** @} */ // end of IsectContainer group
#endif
