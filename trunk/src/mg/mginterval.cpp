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
// MGInterval�N���X�̎���
//

//
//Constructor:�R���X�g���N�^
// void �R���X�g���N�^
MGInterval::MGInterval():m_high(0.0),m_low(1.+MGTolerance::wc_zero()){;}

// Interval�^�C�v�ƈ�_���w�肵��Interval�𐶐��B
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
//����l��ԋp����B
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
//  �C���^�[�o�����L���̎��L���ŉ����l��ԋp����B
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

//  ���g��Interval�̏���l��ύX����B
void MGInterval::set_high_point(const MGEReal& high_point)
{	m_high = high_point;}

//  ���g��Interval�̉����l��ύX����B
void MGInterval::set_low_point(const MGEReal& low_point)
{    m_low = low_point;}

// �Q��
//  Interval�� empty ���ǂ������肷��B
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

//  Interval���L�����ǂ������肷��B
bool MGInterval::finite() const{
	if(empty()) return true;
    return m_low.finite()&&m_high.finite();
}

//  Interval������L�����ǂ������肷��B
bool MGInterval::finite_above() const{
    return m_high.finite();
}
 
//  Interval�������L�����ǂ������肷��B
bool MGInterval::finite_below () const{
    return m_low.finite();
}

//  Interval���������ǂ������肷��B
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

//  Interval�̕�Ԃ�ԋp����B�^����ꂽ param �ɑ΂��āA
//  (1-param)*low_point()+param*high_point() ��ԋp����B
//  ( 0 > param, 1 < param �̎��͊O�� )
double MGInterval::interpolate (double param) const{
	return (1.0-param)*m_low.value()+param*m_high.value();
}

//  ��[�Ɖ��[�̍��ق�ԋp����Bempty�̂Ƃ��̓[���܂��͕��̒l��ԋp����B
MGEReal MGInterval::length() const{
      return (m_high-m_low);
}   

//  Interval���L���̎��L���Œ��_��ԋp����B
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
//  �C���^�[�o���̃^�C�v��ԋp����B
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
// ���Z�q�̑��d��`
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

//  ���g��Interval�Ɨ^����ꂽInterval�̉��Z���s��object
//  �𐶐�����BSame as operator|.
MGInterval MGInterval::operator+ (const MGInterval& intv1) const{
     return (*this)|intv1;
}

//  ���g��Interval�� double �̉��Z���s���u�W�F�N�g��������B
MGInterval MGInterval::operator+ (double value) const{
     MGInterval intv1 = *this ;
	 intv1 += value;
     return intv1;
}

//  ���g��Interval�ɉ��Z�����g��Interval�Ƃ���B
MGInterval& MGInterval::operator+= (const MGInterval& intv1){
	*this |= intv1;
    return *this;
}

MGInterval& MGInterval::operator+= (double value){
	m_low+=value; m_high+=value;
    return *this;
}

//  �O�u�P���}�C�i�X�Bobject�𐶐��B
MGInterval MGInterval::operator- () const{
    MGInterval intv1=*this;
	intv1.m_low=-m_high; intv1.m_high=-m_low;
	return intv1;
}  

//  ���g��Interval�Ɨ^����ꂽInterval�̌��Z���s��object
//  �𐶐�����B
MGInterval MGInterval::operator- (const MGInterval& intv1) const{
    MGInterval intv2 = *this;
    intv2 -=intv1;
    return intv2;
}

//  ���g��Interval�� double �̌��Z���s��Object�𐶐�����B
MGInterval MGInterval::operator- (double value) const{
    MGInterval intv2 = *this;
    intv2 += (-value);
    return intv2;
}

//  ���g��Interval�����Z�����g��Interval�Ƃ���B
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

//  �X�J���[�̏�Z���s��object�𐶐�����B
MGInterval MGInterval::operator* (double value) const{
    MGInterval intv1 = *this;
	intv1 *= value;
    return intv1;
}         

//  �X�J���[�̏�Z���s�����g��Interval�Ƃ���B
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

//  �X�J���[�̏��Z���s��object�𐶐�����B
MGInterval MGInterval::operator/ (double value) const{
    MGInterval intv1 = *this;
	intv1 /= value;
    return intv1;
}

//  �X�J���[�̏��Z���s�����g��Interval�Ƃ���B
MGInterval & MGInterval::operator/= (double value){
	double a=1./value;
    return (*this *= a);
}

//  ���g��Interval�Ɨ^����ꂽInterval����������Interval��
//  ��������B
MGInterval MGInterval::operator| (const MGInterval& i2) const{
     MGInterval i1 = *this;
	 i1 |= i2;
     return i1;
}

//  ���g��Interval�Ɨ^����ꂽInterval���������Ď��g��Interval
//  �Ƃ���B
MGInterval& MGInterval::operator|= (const MGInterval& i2){
	if(empty()){ m_low=i2.m_low; m_high=i2.m_high;}
	else{
		if(m_high<i2.m_high) m_high=i2.m_high;
		if(m_low>i2.m_low) m_low=i2.m_low;
	}
     return *this;
}         

//  ���g��Interval�Ɨ^����ꂽInterval�̋��ʕ�����Interval��
//  ��������B
MGInterval MGInterval::operator& (const MGInterval& i2) const{
     MGInterval i1 = *this;
     i1 &= i2;
     return i1;
}

//  ���g��Interval�Ɨ^����ꂽInterval���ʕ��������g��Interval
//  �Ƃ���B
MGInterval& MGInterval::operator&= (const MGInterval& i2){
	if(i2.m_low>m_low) m_low=i2.m_low;
	if(i2.m_high<m_high) m_high=i2.m_high;
    return *this;
}

//  Boolean ���Z

//  ���g��Interval�Ɨ^����ꂽInterval�����ʕ����������Ă��邩
//  �ǂ�����ԋp����B
bool MGInterval::operator&& (const MGInterval& i2)const{
    MGInterval temp = *this & i2;
    return !(temp.empty());
}

//  ���g��Interval�� empty �̏ꍇ False(0) ��ԋp���A�^����ꂽInterval
//  �� empty �̏ꍇ True(1) ��ԋp����B�^����ꂽInterval�̉��[�����g
//  ��Interval�̉��[���傫���A�^����ꂽInterval�̏�[�����g��
//  Interval�̏�[������������ True(1) ��ԋp����B
//  �܂��A������Interval�� empty �̎� True(1) ��ԋp����B
bool MGInterval::operator>> (const MGInterval& i2) const{
	if(i2.empty()) return 1;
	else if(empty()) return 0;
	else return (m_low<=i2.m_low && m_high>=i2.m_high);
}
 
//  �^����ꂽ�l�����g��Interval���ɂ��邩�ǂ����ԋp����B�܂܂��
//  �� True(1) ��ԋp����B
bool MGInterval::operator>> (const MGEReal& value) const{
	return (value>=m_low && value<=m_high);
}

//  ���g��Interval�� empty �̏ꍇ True(1) ��ԋp���A�^����ꂽInterval
//  �� empty �̏ꍇ False(0) ��ԋp����B�^����ꂽInterval�̉��[�����g
//  ��Interval�̉��[��菬�����A�^����ꂽInterval�̏�[�����g��
//  Interval�̏�[�����傫���� True(1) ��ԋp����B
//  �܂��A������Interval�� empty �̎� False(0) ��ԋp����B
bool MGInterval::operator<< (const MGInterval& i2) const{
	return i2>>(*this);
}

//  �^����ꂽ�l�����g��Interval�͈͊O�ɂ��邩�ǂ����ԋp����B�͈͊O��
//  �� True(1) ��ԋp����B���g��Interval�� empty �̏ꍇ�� True(1) ��
//  �ԋp����B
bool MGInterval::operator<< (const MGEReal& value) const{
	return !((*this)>>value);
}

//  �Q��Interval�����ꂩ�ǂ����𔻒肷��B
//  ����ł��鎞 True(1) ��ԋp
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

//  �Q��Interval�����ꂩ�ǂ����𔻒肷��B
//  ����łȂ� �� True(1) ��ԋp
bool MGInterval::operator != (const MGInterval& i2) const{
	return !(*this==i2);
}

//  �Z�p�I��r
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

//  ���g��Interval�� double �̉��Z���s���u�W�F�N�g��������B
MGInterval operator +(double data, const MGInterval & i2){
     MGInterval i1 = i2;
     i1 += data;
     return i1;
}

//  �X�J���[�̏�Z���s��object�𐶐�����B
MGInterval operator *(double data, const MGInterval &i2){
	MGInterval i1 = i2;
	i1 *= data;
	return i1;
}
