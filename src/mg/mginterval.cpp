/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Interval.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGInterval.cc
// MGIntervalクラスの実装
//

//
//Constructor:コンストラクタ
// void コンストラクタ
MGInterval::MGInterval():m_high(0.0),m_low(1.+MGTolerance::wc_zero()){;}

// Intervalタイプと一点を指定してIntervalを生成。
MGInterval::MGInterval(
	MGINTERVAL_TYPE type,
	double point        ){
	 switch ( type ) {
	   case MGINTERVAL_FINITE:
            m_high=m_low=MGEReal(point);
            break;
 	   case MGINTERVAL_FINITE_ABOVE:
            m_high=MGEReal(point);m_low=MGEReal(MGINFINITE_MINUS);
            break;
	   case MGINTERVAL_FINITE_BELOW:
            m_high=MGEReal(MGINFINITE_PLUS);m_low=MGEReal(point);
            break;
	   case MGINTERVAL_INFINITE :
            m_high=MGEReal(MGINFINITE_PLUS);m_low=MGEReal(MGINFINITE_MINUS);
            break;
	   case MGINTERVAL_EMPTY : 
            m_high=MGEReal(0.0);m_low=MGEReal(1.+MGTolerance::wc_zero());
	 }
}

//Construct interval from two MGEReal points.
MGInterval::MGInterval(const MGEReal& point1,
					   const MGEReal& point2)
:m_low(point1),m_high(point2){;}

// Construct an interval which contains both input interval and a point.
MGInterval::MGInterval(const MGInterval& intvl,const MGEReal& p)
:m_low(intvl.m_low), m_high(intvl.m_high){
	if(empty()){ m_low=p; m_high=p;}
	else{ if(p<m_low) m_low=p; if(p>m_high) m_high=p;}
}

//Construct a interval of one point.
MGInterval::MGInterval(double t):m_low(t), m_high(t){;}

//
// Member Function
//

//Chnage the interval range to [t0, t1] if t0<t1,
//                          to [t1, t0] if t1<t0.
void MGInterval::change_range(
	double t0, double t1	//new range.
){
	if(t0<t1){
		set_low_point(t0); set_high_point(t1);
	}else{
		set_low_point(t1); set_high_point(t0);
	}
}

//Return high point of the interval.
//上限値を返却する。
double MGInterval::high_point() const
{if(empty()) return 0.0; else return m_high.value();}

//Test if this finite interval includes the value, taking account relative
//error of the interval into the account.
//This must be an finite interval.
bool MGInterval::includes(double t)const{
	double low=m_low.value(), high=m_high.value();
	double len=high-low;
	double error=len*MGTolerance::rc_zero();
	if(t<low-error) return false;
	if(t>high+error) return false;
	return true;
}

//Return low point of the interval.
//  インターバルが有限の時有効で下限値を返却する。
double MGInterval::low_point() const
{if(empty()) return 0.0; else return m_low.value();}

//Compute relative error of the interval.
double MGInterval::relative_error() const{
	MGEReal span=m_high-m_low;
	if(span.finite()) return span.value()*MGTolerance::rc_zero();
    else              return MGTolerance::wc_zero();
}

//Set this as null(empty) interval
void MGInterval::set_null(double t){
	set_low_point(t);
	set_high_point(t-MGTolerance::rc_zero());
}

//  自身のIntervalの上限値を変更する。
void MGInterval::set_high_point(const MGEReal& high_point)
{	m_high = high_point;}

//  自身のIntervalの下限値を変更する。
void MGInterval::set_low_point(const MGEReal& low_point)
{    m_low = low_point;}

// 参照
//  Intervalが empty かどうか判定する。
//***Interval of m_low==m_high(including error) is not empty.***
bool MGInterval::empty() const{
	if(m_high.minus_infinite()) return true;
	else if(m_low.plus_infinite()) return true;
    else return ((m_low-MGTolerance::wc_zero())>m_high);
}

//Expand the interval so that this includes the doube val.
void MGInterval::expand(const MGEReal& val){
	if(empty()){ m_low=m_high=val;}
	else{ if(val<m_low) m_low=val; if(val>m_high) m_high=val;}
}

//  Intervalが有限かどうか判定する。
bool MGInterval::finite() const{
	if(empty()) return true;
    return m_low.finite()&&m_high.finite();
}

//  Intervalが上方有限かどうか判定する。
bool MGInterval::finite_above() const{
    return m_high.finite();
}
 
//  Intervalが下方有限かどうか判定する。
bool MGInterval::finite_below () const{
    return m_low.finite();
}

//  Intervalが無限かどうか判定する。
bool MGInterval::infinite() const{
    return m_low.minus_infinite() && m_high.plus_infinite();
}

//Test if infinite_above.
bool MGInterval::infinite_above() const{
    return m_high.plus_infinite();
}
//Test if infinite_below.
bool MGInterval::infinite_below() const{
    return m_low.minus_infinite();
}

//  Intervalの補間を返却する。与えられた param に対して、
//  (1-param)*low_point()+param*high_point() を返却する。
//  ( 0 > param, 1 < param の時は外分 )
double MGInterval::interpolate (double param) const{
	return (1.0-param)*m_low.value()+param*m_high.value();
}

//  上端と下端の差異を返却する。emptyのときはゼロまたは負の値を返却する。
MGEReal MGInterval::length() const{
      return (m_high-m_low);
}   

//  Intervalが有限の時有効で中点を返却する。
double MGInterval::mid_point() const{
	return (m_low.value()+m_high.value())/2.;
}   

//Round the input value into this interval's range, that is:
//round_into_interval(t)=(*this)[0] if t<=(*this)[0],
//round_into_interval(t)=(*this)[1] if t>=(*this)[1].
double MGInterval::round_into_interval(double t)const{
	if(t<(*this)[0]) return low_point();
	else if(t>(*this)[1]) return high_point();
	return t;
}

//Return interval type.
//  インターバルのタイプを返却する。
MGINTERVAL_TYPE MGInterval::type() const
{
	if(empty())                    return MGINTERVAL_EMPTY;
	if(m_low.minus_infinite()){
		if(m_high.plus_infinite()) return MGINTERVAL_INFINITE;
		else                       return MGINTERVAL_FINITE_ABOVE;
	}else{
		if(m_high.plus_infinite()) return MGINTERVAL_FINITE_BELOW;
		else                       return MGINTERVAL_FINITE;
	}
}

//
// 演算子の多重定義
//

//Return low or high point.
//0<=i<=1, and for i=0, return low point, for i=1, return high point.
MGEReal& MGInterval::operator() (size_t i){
	assert(i<=1);
	if(i) return m_high;
	return m_low;
}
const MGEReal& MGInterval::operator() (size_t i)const{
	assert(i<=1);
	if(i) return m_high;
	return m_low;
}
const MGEReal& MGInterval::operator[] (size_t i)const{
	assert(i<=1);
	if(i) return m_high;
	return m_low;
}
MGEReal& MGInterval::operator[] (size_t i){
	assert(i<=1);
	if(i) return m_high;
	return m_low;
}

//  自身のIntervalと与えられたIntervalの加算を行いobject
//  を生成する。Same as operator|.
MGInterval MGInterval::operator+ (const MGInterval& intv1) const{
     return (*this)|intv1;
}

//  自身のIntervalと double の加算を行いブジェクト生成する。
MGInterval MGInterval::operator+ (double value) const{
     MGInterval intv1 = *this ;
	 intv1 += value;
     return intv1;
}

//  自身のIntervalに加算し自身のIntervalとする。
MGInterval& MGInterval::operator+= (const MGInterval& intv1){
	*this |= intv1;
    return *this;
}

MGInterval& MGInterval::operator+= (double value){
	m_low+=value; m_high+=value;
    return *this;
}

//  前置単項マイナス。objectを生成。
MGInterval MGInterval::operator- () const{
    MGInterval intv1=*this;
	intv1.m_low=-m_high; intv1.m_high=-m_low;
	return intv1;
}  

//  自身のIntervalと与えられたIntervalの減算を行いobject
//  を生成する。
MGInterval MGInterval::operator- (const MGInterval& intv1) const{
    MGInterval intv2 = *this;
    intv2 -=intv1;
    return intv2;
}

//  自身のIntervalと double の減算を行いObjectを生成する。
MGInterval MGInterval::operator- (double value) const{
    MGInterval intv2 = *this;
    intv2 += (-value);
    return intv2;
}

//  自身のIntervalを減算し自身のIntervalとする。
MGInterval& MGInterval::operator-= (const MGInterval& i2){
	if(i2.m_high>m_low){
		if(i2.m_high<=m_high){
			if(i2.m_low<=m_low) m_low=i2.m_high;
			else if(i2.m_high==m_high) m_high=i2.m_low;
		}else if(i2.m_low<=m_low){
			m_low=i2.m_high; m_high=i2.m_low;
		}else if(i2.m_low<=m_high){
			m_high=i2.m_low;
		}
	}
    return *this;
}

MGInterval& MGInterval::operator-= (double value){
    *this += ( - value );
    return *this;
}

//  スカラーの乗算を行いobjectを生成する。
MGInterval MGInterval::operator* (double value) const{
    MGInterval intv1 = *this;
	intv1 *= value;
    return intv1;
}         

//  スカラーの乗算を行い自身のIntervalとする。
MGInterval& MGInterval::operator*= (double value)
{
	if(finite()){
		double low=m_low.value(), high=m_high.value();
		double midpoint=(low+high)*0.5;
		double half=(high-low)*value*0.5;
		if(value>=0.){
			m_low=midpoint-half;
			m_high=midpoint+half;
		}else{
			m_low=-midpoint+half;
			m_high=-midpoint-half;
		}
	}else if(value<0.){
		MGEReal temp=m_low; m_low=m_high*-1.; m_high=temp*-1.;
	}
    return *this;
}         

//  スカラーの除算を行いobjectを生成する。
MGInterval MGInterval::operator/ (double value) const{
    MGInterval intv1 = *this;
	intv1 /= value;
    return intv1;
}

//  スカラーの除算を行い自身のIntervalとする。
MGInterval & MGInterval::operator/= (double value){
	double a=1./value;
    return (*this *= a);
}

//  自身のIntervalと与えられたIntervalを結合したIntervalを
//  生成する。
MGInterval MGInterval::operator| (const MGInterval& i2) const{
     MGInterval i1 = *this;
	 i1 |= i2;
     return i1;
}

//  自身のIntervalと与えられたIntervalを結合して自身のInterval
//  とする。
MGInterval& MGInterval::operator|= (const MGInterval& i2){
	if(empty()){ m_low=i2.m_low; m_high=i2.m_high;}
	else{
		if(m_high<i2.m_high) m_high=i2.m_high;
		if(m_low>i2.m_low) m_low=i2.m_low;
	}
     return *this;
}         

//  自身のIntervalと与えられたIntervalの共通部分のIntervalを
//  生成する。
MGInterval MGInterval::operator& (const MGInterval& i2) const{
     MGInterval i1 = *this;
     i1 &= i2;
     return i1;
}

//  自身のIntervalと与えられたInterval共通部分を自身のInterval
//  とする。
MGInterval& MGInterval::operator&= (const MGInterval& i2){
	if(i2.m_low>m_low) m_low=i2.m_low;
	if(i2.m_high<m_high) m_high=i2.m_high;
    return *this;
}

//  Boolean 演算

//  自身のIntervalと与えられたIntervalが共通部分を持っているか
//  どうかを返却する。
bool MGInterval::operator&& (const MGInterval& i2)const{
    MGInterval temp = *this & i2;
    return !(temp.empty());
}

//  自身のIntervalが empty の場合 False(0) を返却し、与えられたInterval
//  が empty の場合 True(1) を返却する。与えられたIntervalの下端が自身
//  のIntervalの下端より大きく、与えられたIntervalの上端が自身の
//  Intervalの上端よりも小さい時 True(1) を返却する。
//  また、両方のIntervalが empty の時 True(1) を返却する。
bool MGInterval::operator>> (const MGInterval& i2) const{
	if(i2.empty()) return 1;
	else if(empty()) return 0;
	else return (m_low<=i2.m_low && m_high>=i2.m_high);
}
 
//  与えられた値が自身のInterval内にあるかどうか返却する。含まれる
//  時 True(1) を返却する。
bool MGInterval::operator>> (const MGEReal& value) const{
	return (value>=m_low && value<=m_high);
}

//  自身のIntervalが empty の場合 True(1) を返却し、与えられたInterval
//  が empty の場合 False(0) を返却する。与えられたIntervalの下端が自身
//  のIntervalの下端より小さく、与えられたIntervalの上端が自身の
//  Intervalの上端よりも大きい時 True(1) を返却する。
//  また、両方のIntervalが empty の時 False(0) を返却する。
bool MGInterval::operator<< (const MGInterval& i2) const{
	return i2>>(*this);
}

//  与えられた値が自身のInterval範囲外にあるかどうか返却する。範囲外の
//  時 True(1) を返却する。自身のIntervalが empty の場合も True(1) を
//  返却する。
bool MGInterval::operator<< (const MGEReal& value) const{
	return !((*this)>>value);
}

//  ２つのIntervalが同一かどうかを判定する。
//  同一である時 True(1) を返却
bool MGInterval::operator== (const MGInterval& i2) const{
	if(finite() && i2.finite()){
		double mzero=MGTolerance::mach_zero();
		double low1=low_point(), high1=high_point();
		double low2=i2.low_point(), high2=i2.high_point();
		double span=high1-low1;
		double span2=high2-low2;
		if(span<span2) span=span2;
		if(span<=mzero) return MGAEqual(low1,low2);
		return MGRZero2(low1-low2,span)&&MGRZero2(high1-high2,span);
	}else return (m_low==i2.m_low && m_high==i2.m_high);
}

//  ２つのIntervalが同一かどうかを判定する。
//  同一でない 時 True(1) を返却
bool MGInterval::operator != (const MGInterval& i2) const{
	return !(*this==i2);
}

//  算術的比較
bool MGInterval::operator < (const MGInterval& i2) const{
	return (m_low<i2.m_low && m_high<i2.m_high);
}
bool MGInterval::operator< (const MGEReal& value) const{
	return (m_high<value);
}
bool MGInterval::operator> (const MGInterval& i2) const{
	return (m_low>i2.m_low && m_high>i2.m_high);
}
bool MGInterval::operator> (const MGEReal& value) const{
	return (m_low>value);
}
bool MGInterval::operator<= (const MGInterval& i2) const{
	return *this==i2 || *this<i2;
}
bool MGInterval::operator<= (const MGEReal& value) const{
	double error=relative_error();
	double save=MGTolerance::set_wc_zero(error);
	bool result=(m_high<=value);
	MGTolerance::set_wc_zero(save);
	return result;
}
bool MGInterval::operator>= (const MGInterval& i2) const{
	return *this==i2 || *this>i2;
}         
bool MGInterval::operator>= (const MGEReal& value) const{
	double save=MGTolerance::set_wc_zero(relative_error());
	bool result=(m_low>=value);
	MGTolerance::set_wc_zero(save);
	return result;
}

bool operator>(const MGEReal& t, const MGInterval& i){return i<t;}
bool operator<(const MGEReal& t, const MGInterval& i){return i>t;}
bool operator>=(const MGEReal& t, const MGInterval& i){return i<=t;}
bool operator<=(const MGEReal& t, const MGInterval& i){return i>=t;}

bool MGInterval::operator< (double t) const{
	if(m_high.plus_infinite()) return false;
	return m_high.m_value<t;
}
bool operator>(double t, const MGInterval& i){return i<t;}
bool MGInterval::operator> (double t) const{
	if(m_low.minus_infinite()) return false;
	return m_low.m_value>t;
}
bool operator<(double t, const MGInterval& i){return i>t;}

//Global Function

//  自身のIntervalと double の加算を行いブジェクト生成する。
MGInterval operator +(double data, const MGInterval & i2){
     MGInterval i1 = i2;
     i1 += data;
     return i1;
}

//  スカラーの乗算を行いobjectを生成する。
MGInterval operator *(double data, const MGInterval &i2){
	MGInterval i1 = i2;
	i1 *= data;
	return i1;
}
