/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGInterval_HH_
#define _MGInterval_HH_
#include "mg/EReal.h"
/** @addtogroup BASE
 *  @{
 */

//  MGInterval.h
//  header for class MGInterval
// Forward Declaration
class MGIfstream;
class MGOfstream;

/// Interval of 1 dimension, i.e. MGInterval is a real line.
///  数直線上の区間を表す。２つの MGEReal で表現される。
class MGCLASS MGInterval{

public:

MGDECL friend bool operator>(const MGEReal&, const MGInterval&);
MGDECL friend MGInterval operator+ (double, const MGInterval& );

///Scalar multiplication to the interval.
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
MGDECL friend MGInterval operator* (double a, const MGInterval& );

MGDECL friend bool operator<(const MGEReal&, const MGInterval&);
MGDECL friend bool operator>=(const MGEReal&, const MGInterval&);
MGDECL friend bool operator<=(const MGEReal&, const MGInterval&);
MGDECL friend bool operator>(double, const MGInterval& i);
MGDECL friend bool operator<(double, const MGInterval& i);

///String stream function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGInterval&);

//////////////Constructor////////////

/// void Constructor
MGInterval();

///Construct interval from interval type and a point.
///When the type is MGINTERVAL_FINITE, start and end is the point.
///MGINTERVAL_TYPE:MGINTERVAL_EMPTY, MGINTERVAL_FINITE,
///MGINTERVAL_FINITE_ABOVE, MGINTERVAL_FINITE_BELOW, or MGINTERVAL_INFINITE.
///MGINTERVAL_FINITE_ABOVE is infinite below and finite above.
///MGINTERVAL_FINITE_BELOW is infinite above and finite below.
/// インターバルタイプと一点を指定してインターバルを生成。
///     MGINTERVAL_FINITE 指定時 指定された一点の点範囲のインターバル
///     MGINTERVAL_FINITE_ABOVE, MGINTERVAL_FINITE_BELOW 指定時 
///     有限の一点を指定
///     MGINTERVAL_INFINITE 指定時 double の値の指定は必要なし
explicit MGInterval(MGINTERVAL_TYPE, double=0.0);

///Construct interval from two MGEReal points.
///t1 is lower limit, and t2 is upper limit.
///When t1> t2, the interval is empry.
/// ２つのMGERealの値からインターバルを生成。
///***** This is the fundamental constructor.
MGInterval(const MGEReal& t1, const MGEReal& t2); 

/// Construct an interval which contains both input interval and a value.
MGInterval(const MGInterval&, const MGEReal&); 

///Construct a interval of one point.
explicit MGInterval(double t);

////////////Operator overload////////////

///Return low or high point.
///0<=i<=1, and for i=0, return low point, for i=1, return high point.
const MGEReal& operator() (size_t i)const;
MGEReal& operator() (size_t i);
const MGEReal& operator[] (size_t i)const;
MGEReal& operator[] (size_t i);

///Addition of two intervals. Generated interval is the minimum one
///that includes both intervals. The same as operator| .
///  自身のインターバルと与えられたインターバルの加算を行いオブジェクト
///  を生成する。
MGInterval operator+ (const MGInterval& ) const;

///Translation of the interval.
///  自身のインターバルと double の加算を行いブジェクト生成する。
MGInterval operator+ (double) const;

///Addition of two intervals. Updated interval is the minimum one
///that includes both original intervals. The same as operator|=. 
///  自身のインターバルに加算し自身のインターバルとする。
MGInterval& operator+= (const MGInterval& );

///Translation of the interval.
MGInterval& operator+= (double);

///Unary minus.
///  前置単項マイナス。オブジェクトを生成。
MGInterval operator- () const;

///Subtraction. Subtraction is defined as:
/// "this" minus the common of the two intervals. 
///  自身のインターバルと与えられたインターバルの減算を行いオブジェクト
///  を生成する。
MGInterval operator- (const MGInterval& ) const ;

///Translation of the interval.
///  自身のインターバルと double の減算を行いブジェクト生成する。
MGInterval operator- (double) const;

///Subtraction. Subtraction is defined as:
/// "this" minus the common of the two intervals. 
///  自身のインターバルに減算し自身のインターバルとする。
MGInterval& operator-= (const MGInterval& );

///Translation of the interval.
MGInterval& operator-= (double);

///Scalar multiplication to the interval.
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
///  スカラーの乗算を行いオブジェクトを生成する。
MGInterval operator* (double a) const;

///Scalar multiplication to the interval.
///  スカラーの乗算を行い自身のインターバルとする。
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
MGInterval& operator*= (double a);

///Scalar division to the interval.
///  スカラーの除算を行いオブジェクトを生成する。
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the division of the value a.
MGInterval operator/ (double a) const; 

///Scalar division to the interval.
///  スカラーの除算を行い自身のインターバルとする。
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the division of the value a.
MGInterval& operator/= (double a);

///Addition of two intervals. Generated interval is the minimum one
///that includes both intervals. The same as operator+ .
///  自身のインターバルと与えられたインターバルをOR結合したインターバルを
///  生成する。
MGInterval operator| (const MGInterval&) const ;

///Addition of two intervals. Updated interval is the minimum one
///that includes both original intervals. The same as operator+=. 
///自身のインターバルと与えられたインターバルをOR結合して自身のインターバル
///とする。
MGInterval& operator|= (const MGInterval&);

///And operation. Generate the interval common to the both intervals.
///自身のインターバルと与えられたインターバルの共通部分のインターバルを
///生成する。
MGInterval operator& (const MGInterval&) const ;

///And operation. Update to the interval common to the both intervals.
///  自身のインターバルと与えられたインターバル共通部分を自身のインターバル
///  とする。
MGInterval& operator&= (const MGInterval&);

///  Boolean 演算

///Test if the two intervals have common part.
///  自身のインターバルと与えられたインターバルが共通部分を持っているか
///  どうかを返却する。
bool operator&& (const MGInterval&)const;

///Test if this includes intrval2. 
///  自身のインターバルが empty の場合 false(0) を返却し、与えられたインターバル
///  が empty の場合 true(1) を返却する。与えられたインターバルの下端が自身
///  のインターバルの下端より大きく、与えられたインターバルの上端が自身の
///  インターバルの上端よりも小さい時 true(1) を返却する。
///  また、両方のインターバルが empty の時 true(1) を返却する。
bool operator>> (const MGInterval& intrval2) const;

///Test if this includes the value t.
///  与えられた値が自身のインターバル内にあるかどうか返却する。含まれる
///  時 true(1) を返却する。自身のインターバルが empty の場合は false(0) を
///  返却する。
bool operator>> (const MGEReal& t) const;

///Test if this is included in intrval2.
///  自身のインターバルが empty の場合 true(1) を返却し、与えられたインターバル
///  が empty の場合 false(0) を返却する。与えられたインターバルの下端が自身
///  のインターバルの下端より小さく、与えられたインターバルの上端が自身の
///  インターバルの上端よりも大きい時 true(1) を返却する。
///  また、両方のインターバルが empty の時 false(0) を返却する。
bool operator<< (const MGInterval& intrval2) const;

///Test if this excludes the value t.
///  与えられた値が自身のインターバル範囲外にあるかどうか返却する。範囲外の
///  時 true(1) を返却する。自身のインターバルが empty の場合も true(1) を
///  返却する。
bool operator<< (const MGEReal& t) const;

///Test if two intervals are equal(tolerance included)
///  ２つのインターバルが同一かどうかを判定する。
///  同一である時 true(1) を返却
bool operator== (const MGInterval&) const;

///Test if two intervals are not equal(tolerance included)
///  ２つのインターバルが同一かどうかを判定する。
///  同一でない 時 true(1) を返却
bool operator!= (const MGInterval&) const;

///Arithmatic comparison operation:算術的比較
bool operator< (const MGInterval&) const;
bool operator> (const MGInterval&) const;
bool operator<= (const MGInterval&) const;
bool operator>= (const MGInterval&) const;

bool operator< (const MGEReal& ) const;
bool operator> (const MGEReal&) const;
bool operator<= (const MGEReal&) const;
bool operator>= (const MGEReal&) const;

bool operator< (double t) const;
bool operator> (double t) const;

////////////Member Function////////////

///Chnage the interval range to [t0, t1] if t0<t1,
///                          to [t1, t0] if t1<t0.
void change_range(
	double t0, double t1	///<new range.
);

/// Reference
///Test if empty:インターバルが empty かどうか判定する。
///***Interval of m_low==m_high(including error) is not empty.***
bool empty() const;

///Expand the interval so that this includes the doube val.
void expand(const MGEReal& val);

///Test if finite.
///  インターバルが有限かどうか判定する。
bool finite() const;

///Test if finite above(infinite below).
///  インターバルが上方有限かどうか判定する。
bool finite_above() const;
 
///Test if finite below(infinite above).
///  インターバルが下方有限かどうか判定する。
bool finite_below() const;

const MGEReal& high() const {return m_high;};
MGEReal& high() {return m_high;};

///Return high point of the interval. When infinite_above(),
///mgInfiniteVal will be returned.　When empty, return 0.
///上限値を返却する。
double high_point() const;

///Test if this finite interval includes the value, taking account relative
///error of the interval into the account.
///This must be an finite interval.
bool includes(double t)const;

///Test if infinite.
///  インターバルが無限かどうか判定する。
bool infinite() const;

///Test if infinite_above.
bool infinite_above() const;

///Test if infinite_below.
bool infinite_below() const;

///Compute interpolated point of the interval. Given parameter t,
///the interpolated point=(1-t)*low_point()+t*high_point().
///t may be t<0 or t>1. Valid only when finite().
///  インターバルの補間を返却する。与えられた param に対して、
///  (1-t)*low_point()+t*high_point() を返却する。
///  ( 0 > t, 1 < t の時は外分 )
double interpolate(double t) const;

///Test if this is empty interval.
bool is_null() const{return empty();};

///Compute the length of the interval.
///When empty, length()<=0.0.
///上端と下端の差異を返却する。empty は非正値（ゼロまたは負の値）
///を返却する。
MGEReal length() const;

const MGEReal& low() const {return m_low;};
MGEReal& low() {return m_low;};

///Return low point of the interval. When infinite_below(),
/// -mgInfiniteVal will be returned. When empty, return 0.
///  インターバルが有限の時有効で下限値を返却する。空の時 0.0 を
///  返却する。
double low_point() const;

///Return mid point of the interval. Valid only when finite().
///  インターバルが有限の時有効で中点を返却する。
double mid_point() const;

///Compute relative_error*length();
double relative_error() const;

///Round the input value into this interval's range, that is:
///round_into_interval(t)=(*this)[0] if t<=(*this)[0],
///round_into_interval(t)=(*this)[1] if t>=(*this)[1].
double round_into_interval(double t)const;

///Set this as null(empty) interval
void set_null(double t=0.);
void set_empty(double t=0.){set_null(t);};

///Update high point of the interval. When given point is less than
///low point, the interval becomes empty.
///  自身のインターバルの上限値を変更する。与えられた値が下限値より小さい
///  場合、自身のインターバルは empty になり、上下限値の入れ換えはない。
void set_high_point(const MGEReal&);
void set_high(const MGEReal& t){set_high_point(t);};

///Update low point of the interval. When given point is greater than
///high point, the interval becomes empty.
///  自身のインターバルの下限値を変更する。与えられた値が上限値より大きい
///  場合、自身のインターバルは empty になり、上下限値の入れ換えはない。
void set_low_point(const MGEReal&);
void set_low(const MGEReal& t){set_low_point(t);};

///Return interval type.
///  インターバルのタイプを返却する。
MGINTERVAL_TYPE type() const;

///Dump Functions.
///Calculate dump size
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////// Member data.////////////

	///  Two points on real line. 数直線上の２点 
	MGEReal m_high;	///<larger one.
	MGEReal m_low;	///<smaller one.

};

/** @} */ // end of BASE group
#endif
