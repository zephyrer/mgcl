/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/EReal.h"
#include "mg/Interval.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGEReal.cpp
//
// MGEReal is extended real number, i.e., it includes minus infinite and 
// plus infinite except ordinary real number.

MGEReal::MGEReal(MGINFINITE_TYPE infinite)	//infinite=-1 means minus_infinite,
											//         +1 means plus_infinite.
:m_value(0.0){
	if(infinite==MGINFINITE_MINUS)
		set_minus_infinite();
	else if(infinite==MGINFINITE_PLUS)
		set_plus_infinite();
}

//Destructor
//	~MGEReal();	

//Member Function

//Operator overload.
MGEReal MGEReal::operator+ (double t) const{
	if(infinite())
		return *this;

	MGEReal er(*this);
	er.m_value+=t;
	return er;
}

MGEReal operator+ (double t, const MGEReal& er){
	return er+t;
}

MGEReal MGEReal::operator+ (const MGEReal& er2) const{
	MGEReal er(*this);
	return er+=er2;
}

MGEReal& MGEReal::operator+= (double t){
	if(infinite())
		return *this;

	m_value+=t;
	return *this;
}

MGEReal& MGEReal::operator+= (const MGEReal& er2){
	int infinite1=infinite_coef();
	int infinite2=er2.infinite_coef();
	if(!infinite1 && !infinite2){
		m_value+=er2.m_value;
	}else{
		int tinf=infinite1+infinite2;
		if(tinf==0)
			set_zero();
		else if(!infinite1)
			m_value=er2.m_value;
	}
	return *this;
}

MGEReal MGEReal::operator-()const{				//Unary minus.
	MGEReal er(*this);
	er.invert();
	return er;
}

MGEReal operator- (double t, const MGEReal& er){
	return -er+t;
}

MGEReal MGEReal::operator-(double t)const{
	MGEReal er(*this);
	return er-=t;
}

MGEReal MGEReal::operator-(const MGEReal& er2) const{
	MGEReal er(*this);
	return er+(-er2);
}

MGEReal& MGEReal::operator-=(double t){
	if(!infinite())
		m_value-=t;
	return *this;
}

MGEReal& MGEReal::operator-=(const MGEReal& er2){
	return (*this)+=(-er2);
}

MGEReal MGEReal::operator* (double t) const{
	MGEReal er(*this);
	return er*=t;
}

MGEReal operator*(double t,const MGEReal& er){
	return er*t;
}

MGEReal MGEReal::operator* (const MGEReal& er2) const{
	MGEReal er(*this);
	return er*=er2;
}

MGEReal& MGEReal::operator*=(double t){
	if(infinite()){
		if(t<0.)
			invert();
		else if(t==0.)
			set_zero();
	}else
		m_value*=t;
	return *this;
}

MGEReal& MGEReal::operator*= (const MGEReal& er2){
	int infinite1=infinite_coef();
	int infinite2=er2.infinite_coef();
	if(!infinite2)
		return (*this)*=er2.m_value;
	if(!infinite1){
		double save=m_value;
		m_value=er2.m_value;
		return (*this)*=save;
	}

	int tinf=infinite1*infinite2;
	if(tinf>0)
		set_plus_infinite();
	else
		set_minus_infinite();

	return *this;
}

MGEReal MGEReal::operator/ (double t) const{
	MGEReal er(*this);
	return er/=t;
}

MGEReal operator/ (double t, const MGEReal& er2){ 
	MGEReal er(t);
	return er/er2;
}

MGEReal MGEReal::operator/ (const MGEReal& er2) const{
	MGEReal er(*this);
	return er/=er2;
}

MGEReal& MGEReal::operator/= (double t){
	if(infinite()){
		if(t<0.)
			invert();
	}else{
		if(MGRZero(t)){
			if(m_value*t>0.)
				set_plus_infinite();
			else
				set_minus_infinite();
		}else
			m_value/=t;
	}
	return *this;
}

MGEReal& MGEReal::operator/= (const MGEReal& er2){
	int infinite1=infinite_coef();
	int infinite2=er2.infinite_coef();
	if(!infinite2)
		return (*this)/=er2.m_value;
	if(!infinite1){
		double save=m_value;
		m_value=er2.m_value;
		return (*this)/=save;
	}

	int tinf=infinite1*infinite2;
	if(tinf>0)
		m_value=1.;
	else
		m_value=-1.;

	return *this;
}

bool MGEReal::operator== (double t) const{
	if(infinite())
		return false;
	return MGAEqual(m_value,t);
}
bool MGEReal::operator==(const MGEReal& er2) const{
	int infinite1=infinite_coef();
	int infinite2=er2.infinite_coef();
	if(infinite1!=infinite2)
		return false;
	if(infinite1)
		return true;
	return MGAEqual(m_value,er2.m_value);
}
bool operator==(double t,const MGEReal& er2){
	return er2==t;
}
bool operator!=(double t,const MGEReal& er2){
	return !(operator==(t,er2));
}

bool MGEReal::operator>(double t) const{
	if(plus_infinite())
		return true;
	else if(minus_infinite())
		return false;
	if(MGAEqual(m_value,t))
		return false;
	return m_value>t;
}
bool MGEReal::operator>(const MGEReal& er2) const{
	int infinite2=er2.infinite_coef();
	if(infinite2){
		int infinite1=infinite_coef();
		if(infinite1)
			return(infinite1>infinite2);
		if(infinite2<0)
			return true;
		return false;
	}else
		return (*this)>er2.m_value;
}
bool operator>(double t,const MGEReal& er2){
	return er2<t;
}

bool MGEReal::operator<(double t) const{
	if(plus_infinite())
		return false;
	else if(minus_infinite())
		return true;
	if(MGAEqual(m_value,t))
		return false;
	return m_value<t;
}
bool operator<(double t,const MGEReal& er2){
	return er2>t;
}

bool MGEReal::operator>=(double t) const{
	if(plus_infinite())
		return true;
	if(minus_infinite())
    	return false;
	if(MGAEqual(m_value,t))
		return true;
	return (m_value>=t);
}

bool MGEReal::operator>=(const MGEReal& er2) const{
	if(plus_infinite())
		return true;
	if(er2.minus_infinite())
		return true;
	if(minus_infinite())
    	return false;
	if(er2.plus_infinite())
		return false;
	if(MGAEqual(m_value,er2.m_value))
		return true;
	return (m_value>=er2.m_value);
}
bool operator>= (double t,const MGEReal& er2){
	return (er2<=t);
}
bool operator<= (double t,const MGEReal& er2){
	return er2>=t;
}
bool MGEReal::operator<= (double t) const{
	if(plus_infinite())
		return false;
	if(minus_infinite())
    	return true;
	if(MGAEqual(m_value,t))
		return true;
	return (m_value<=t);
}

bool MGEReal::equal_base(double t, double base)const{
	if(infinite())
		return false;
	return MGREqual_base(t,m_value,base);
}
bool MGEReal::equal_base(const MGEReal& t, double base)const{
	int infinite1=infinite_coef();
	int infinite2=t.infinite_coef();
	if(infinite1!=infinite2)
		return false;
	if(infinite1)
		return true;

	return MGREqual_base(m_value, t.m_value,base);
}

//return -1 if minus_infinite(), 1 if plus_infinite(), else 0.
int MGEReal::infinite_coef()const{
	if(m_value<=(-mgInfiniteVal))
		return -1;
	if(mgInfiniteVal<=m_value)
		return 1;
	return 0;
}
