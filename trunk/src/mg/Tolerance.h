/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTolerance_HH_
#define _MGTolerance_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"
///  Forward Declarations
class MGIfstream;
class MGOfstream;
class MGEReal;

//  MGTolerance.h

#define MG_MAX_TOL_STACK_SIZE 10

//  Defines tolerance.
//

///MGTolerance is a class to hold various tolerance data.
///Tolerances used are follows:
///(Out of these, wc_zero, line_zero are important. If not familiar, forget 
///rc_zero, mach_zero, m_max_knot_ratio, and angle_zero. Defalut values will be OK.)
///(1) wc_zero(world coordinate zero value):
///Two points distance that should be regarded as coincidence in user's world coordinate.
///Typically, intersection of 2 curves are so computed that the i-points be within
///this tolerance.
///(2) line_zero(Line zero in world coordinates):
///Two curves distance that should be regarded as coincidence.
///line_zero>=wc_zero is recommended.
///(3) rc_zero(normalized relative coordinate zero):
///Two points distance that should be regarded as coincidence in normalized (0.,1.) space.
///This is used when the computation is ecexuted other than world coordinates, like
///in parameter space.
///Generally speaking, when the span of the world interedted is L, L*rc_zero=wc_zero.
///(4) angle_zero: Angle that should be regarded as zero angle in radian.
///Generally speaking, mgDBLPAI*rc_zero=angle_zero.
///(5) m_max_knot_ratio:Maximum ratio of neighboring two spans of a knot vector.
///This is used to stabilize the computation to solve linear equations for B-Spline coefficients. 
///(6) mach_zero(machine zero): the smallest floating value to avoid exception of division zero.
class MGCLASS MGTolerance {
private:
	MGTolerance();
	MGTolerance(const MGTolerance&);
	MGTolerance& operator=(const MGTolerance&);

public:

///String stream function.  �f�o�b�O�֐�
MGDECL friend std::ostream& operator << (std::ostream&, const MGTolerance& );

///�g�������X���l�����ė^����ꂽ double �� 0 ���ǂ������ׂ�B
///Test if data is machine zero.
/// Machine Zero version
MGDECL friend bool MGMZero (double data);

///�g�������X���l�����ė^����ꂽ�Q�� double ����v���邩���ׂ�B
///�������� true(non zero) ��ԋp����-----Absolute Version.
///Test if two double is equal in world coordinate.
MGDECL friend bool MGAEqual (double data1, double data2);

///�g�������X���l�����ė^����ꂽ�l���O�����ׂ�-----Absolute Version.
///Test if data is world coordinate zero.
MGDECL friend bool MGAZero(double data);

///Test if difference of two data is less than MGTolerance::rc_zero(). 
MGDECL friend bool MGREqual (double data1, double data2);

///Test if difference of two data is less than MGTolerance::rc_zero()
///after changing data1 and data2 proportionally for the bigger one to be 1.
MGDECL friend bool MGREqual2(double data1, double data2);

///Test if difference of two data is equal.
///Comparison is:
///test if abs(data1-data2)/base_length is less than MGTolerance::rc_zero().
MGDECL friend bool MGREqual_base(double data1, double data2, double base_length);
MGDECL friend bool MGREqual_base(MGEReal data1, MGEReal data2, const MGEReal& base_length);

///�g�������X���l�����ė^����ꂽ�l���O�����ׂ�1-----Relative Version
///Test if data is less or equal to rc_zero().
MGDECL friend bool MGRZero(double data);

///�g�������X���l�����ė^����ꂽ�l���O�����ׂ�2-----Relative Version
///Test if data is less or equal to rc_zero() compared to base_length.
///Comparison is done after data and base_length are so changed
///that base_length is 1.
///If base_length is zero, MGRZero2 returns always false.
MGDECL friend bool MGRZero2(double data, double base_length);
MGDECL friend bool MGRZero2(double data, const MGEReal& base_length);

///  �p�x�̃R�T�C���l����͂��āC�p�x�����p�����ׂ�B
///�i�g�������X���l������)
///Test if angle is right from the value cos(angle).
MGDECL friend bool MGRight_angle(double cos_data);

///  �g�������X���l�����ė^����ꂽ�p�x(Radian) ���O�����ׂ�B
///Test if data is zero angle taking tolerance into account.
/// If data is small enough, data is almost equal to sin(data).
/// This fact says we can use sine value instead of radian as data input.
MGDECL friend bool MGZero_angle (double data);

//////////// Destructor. ���z�f�X�g���N�^ ////////////
	~MGTolerance();

//////////// Member function.  �����o�֐� //////////

	static MGTolerance& instance();

//Update.  �X�V

///Return currently used  tolerance stack length.
size_t stack_length() const{return m_count;};

///Return maximum tolerance stack size.
size_t max_stack_size() const{return MG_MAX_TOL_STACK_SIZE;};

///Update m_mach_zero. Returned is old mach_zero used so far.
///m_mach_zero��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_mach_zero( double );

///Update m_wc_zero. Returned is old wc_zero used so far.
///m_wc_zero��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_wc_zero( double );

///Update m_rc_zero. Returned is old rc_zero used so far.
/// m_rc_zero��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_rc_zero( double );

///Update m_angle_zero. Returned is old angle_zero used so far.
///m_angle_zero��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_angle_zero( double );

///Update m_line_zero. Returned is old line_zero used so far.
///m_line_zero��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_line_zero( double );

///Update m_max_knot_ratio. Returned is old max_knot_ratio used so far.
///m_max_knot_ratio��ύX����B����܂Ŏg�p����Ă����f�[�^���ԋp�����B
static double set_max_knot_ratio( double );

///Push all the tolerance to stack.  �X�^�b�N�� push ����B
static void push( );

///Pop all the tolerance from stack.  �X�^�b�N�� pop ����B
static void pop( );

///Reference.  �Q��
///Return m_mach_zero. m_mach_zero��ԋp����B
static double mach_zero() {return instance().m_mach_zero;};

///Return m_wc_zero.  m_wc_zero��ԋp����B
static double wc_zero() {return instance().m_wc_zero;};

///Return square of m_wc_zero. m_wc_zero �̂Q���ԋp����
static double wc_zero_sqr(){return instance().m_wc_zero_sqr;};

///Return m_rc_zero. m_rc_zero��ԋp����B
static double rc_zero() {return instance().m_rc_zero;};

///Return square of m_wc_zero. m_wc_zero �̂Q���ԋp����
static double rc_zero_sqr(){return instance().m_rc_zero_sqr;};

///Return m_angle_zero.  m_angle_zero��ԋp����B
static double angle_zero(){return instance().m_angle_zero;};

///Return m_line_zero.  m_line_zero��ԋp����B
static double line_zero(){return instance().m_line_zero;};

///Return m_max_knot_ratio.  m_max_knot_ratio��ԋp����B
static double max_knot_ratio(){return instance().m_max_knot_ratio;};

///Dump Functions
unsigned dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////// Memeber data. �����o�f�[�^ //////////

	double m_mach_zero;	///< �L���l�Ƃ��Ċ���Z���ł���ŏ��̐��l
				///<Machine zero
	double m_wc_zero;	///< �������Ƃ݂Ȃ��Q�_�Ԃ̋���(���[���h���W�j
				///<Two points distance that should be regarded as coincidence
				///<in user's world coordinate.
	double m_wc_zero_sqr;///< m_wc_zero �̂Q��
				///<Square of m_wc_zero.
	double m_rc_zero;	///< �������Ƃ݂Ȃ��Q�_�Ԃ̋���(���΍��W�j
				///<Two points distance that should be regarded as coincidence
				///<in normalized (0.,1.) space.
	double m_rc_zero_sqr;///< m_wc_zero �̂Q��
				///<Square of m_rc_zero.
	double m_angle_zero;///< �[���Ƃ݂Ȃ��p�x
				///<Angle that should be regarded as zero angle.
	double m_line_zero;				///<Two curves distance that should be regarded as coincidence.
				///<�Q�Ȑ��𓙂����Ƃ݂Ȃ��Ȑ��Ԃ̋���
	double m_max_knot_ratio;///<Maximum ratio of neighboring two spans of a knot vector.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) >= m_max_knot_ratio, t(i+1) and
		///<t(i) are regarded same(multiple) knots.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) <= 1/m_max_knot_ratio, t(i) and
		///<t(i-1) are regarded same(multiple) knots.
		///< �ׂ荇��Knot�̔�̍ő�l

///Tolerance Stack.    
    int m_count;				   ///< �X�^�b�N�J�E���^
	double m_mach_zero_stack[MG_MAX_TOL_STACK_SIZE];	///< m_mach_zero �X�^�b�N
	double m_wc_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_wc_zero �X�^�b�N
	double m_rc_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_rc_zero �X�^�b�N
	double m_angle_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_angle_zero �X�^�b�N
	double m_line_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_line_zero �X�^�b�N
	double m_max_knot_ratio_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_max_knot_ratio �X�^�b�N

};

/// Compute radian angle from cosine and sine value.
///Function's return value is angle in radian, from zero to 2PAI.
MGEXTERN double MGAngle(double ca		///<cosine value
					  , double);	///<sine value

/** @} */ // end of BASE group
#endif
