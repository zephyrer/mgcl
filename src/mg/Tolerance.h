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

///String stream function.  デバッグ関数
MGDECL friend std::ostream& operator << (std::ostream&, const MGTolerance& );

///トレランスを考慮して与えられた double が 0 かどうか調べる。
///Test if data is machine zero.
/// Machine Zero version
MGDECL friend bool MGMZero (double data);

///トレランスを考慮して与えられた２つの double が一致するか調べる。
///等しい時 true(non zero) を返却する-----Absolute Version.
///Test if two double is equal in world coordinate.
MGDECL friend bool MGAEqual (double data1, double data2);

///トレランスを考慮して与えられた値が０か調べる-----Absolute Version.
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

///トレランスを考慮して与えられた値が０か調べる1-----Relative Version
///Test if data is less or equal to rc_zero().
MGDECL friend bool MGRZero(double data);

///トレランスを考慮して与えられた値が０か調べる2-----Relative Version
///Test if data is less or equal to rc_zero() compared to base_length.
///Comparison is done after data and base_length are so changed
///that base_length is 1.
///If base_length is zero, MGRZero2 returns always false.
MGDECL friend bool MGRZero2(double data, double base_length);
MGDECL friend bool MGRZero2(double data, const MGEReal& base_length);

///  角度のコサイン値を入力して，角度が直角か調べる。
///（トレランスを考慮する)
///Test if angle is right from the value cos(angle).
MGDECL friend bool MGRight_angle(double cos_data);

///  トレランスを考慮して与えられた角度(Radian) が０か調べる。
///Test if data is zero angle taking tolerance into account.
/// If data is small enough, data is almost equal to sin(data).
/// This fact says we can use sine value instead of radian as data input.
MGDECL friend bool MGZero_angle (double data);

//////////// Destructor. 仮想デストラクタ ////////////
	~MGTolerance();

//////////// Member function.  メンバ関数 //////////

	static MGTolerance& instance();

//Update.  更新

///Return currently used  tolerance stack length.
size_t stack_length() const{return m_count;};

///Return maximum tolerance stack size.
size_t max_stack_size() const{return MG_MAX_TOL_STACK_SIZE;};

///Update m_mach_zero. Returned is old mach_zero used so far.
///m_mach_zeroを変更する。それまで使用されていたデータが返却される。
static double set_mach_zero( double );

///Update m_wc_zero. Returned is old wc_zero used so far.
///m_wc_zeroを変更する。それまで使用されていたデータが返却される。
static double set_wc_zero( double );

///Update m_rc_zero. Returned is old rc_zero used so far.
/// m_rc_zeroを変更する。それまで使用されていたデータが返却される。
static double set_rc_zero( double );

///Update m_angle_zero. Returned is old angle_zero used so far.
///m_angle_zeroを変更する。それまで使用されていたデータが返却される。
static double set_angle_zero( double );

///Update m_line_zero. Returned is old line_zero used so far.
///m_line_zeroを変更する。それまで使用されていたデータが返却される。
static double set_line_zero( double );

///Update m_max_knot_ratio. Returned is old max_knot_ratio used so far.
///m_max_knot_ratioを変更する。それまで使用されていたデータが返却される。
static double set_max_knot_ratio( double );

///Push all the tolerance to stack.  スタックを push する。
static void push( );

///Pop all the tolerance from stack.  スタックを pop する。
static void pop( );

///Reference.  参照
///Return m_mach_zero. m_mach_zeroを返却する。
static double mach_zero() {return instance().m_mach_zero;};

///Return m_wc_zero.  m_wc_zeroを返却する。
static double wc_zero() {return instance().m_wc_zero;};

///Return square of m_wc_zero. m_wc_zero の２乗を返却する
static double wc_zero_sqr(){return instance().m_wc_zero_sqr;};

///Return m_rc_zero. m_rc_zeroを返却する。
static double rc_zero() {return instance().m_rc_zero;};

///Return square of m_wc_zero. m_wc_zero の２乗を返却する
static double rc_zero_sqr(){return instance().m_rc_zero_sqr;};

///Return m_angle_zero.  m_angle_zeroを返却する。
static double angle_zero(){return instance().m_angle_zero;};

///Return m_line_zero.  m_line_zeroを返却する。
static double line_zero(){return instance().m_line_zero;};

///Return m_max_knot_ratio.  m_max_knot_ratioを返却する。
static double max_knot_ratio(){return instance().m_max_knot_ratio;};

///Dump Functions
unsigned dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////// Memeber data. メンバデータ //////////

	double m_mach_zero;	///< 有効値として割り算ができる最小の数値
				///<Machine zero
	double m_wc_zero;	///< 等しいとみなす２点間の距離(ワールド座標）
				///<Two points distance that should be regarded as coincidence
				///<in user's world coordinate.
	double m_wc_zero_sqr;///< m_wc_zero の２乗
				///<Square of m_wc_zero.
	double m_rc_zero;	///< 等しいとみなす２点間の距離(相対座標）
				///<Two points distance that should be regarded as coincidence
				///<in normalized (0.,1.) space.
	double m_rc_zero_sqr;///< m_wc_zero の２乗
				///<Square of m_rc_zero.
	double m_angle_zero;///< ゼロとみなす角度
				///<Angle that should be regarded as zero angle.
	double m_line_zero;				///<Two curves distance that should be regarded as coincidence.
				///<２曲線を等しいとみなす曲線間の距離
	double m_max_knot_ratio;///<Maximum ratio of neighboring two spans of a knot vector.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) >= m_max_knot_ratio, t(i+1) and
		///<t(i) are regarded same(multiple) knots.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) <= 1/m_max_knot_ratio, t(i) and
		///<t(i-1) are regarded same(multiple) knots.
		///< 隣り合うKnotの比の最大値

///Tolerance Stack.    
    int m_count;				   ///< スタックカウンタ
	double m_mach_zero_stack[MG_MAX_TOL_STACK_SIZE];	///< m_mach_zero スタック
	double m_wc_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_wc_zero スタック
	double m_rc_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_rc_zero スタック
	double m_angle_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_angle_zero スタック
	double m_line_zero_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_line_zero スタック
	double m_max_knot_ratio_stack[ MG_MAX_TOL_STACK_SIZE ]; ///< m_max_knot_ratio スタック

};

/// Compute radian angle from cosine and sine value.
///Function's return value is angle in radian, from zero to 2PAI.
MGEXTERN double MGAngle(double ca		///<cosine value
					  , double);	///<sine value

/** @} */ // end of BASE group
#endif
