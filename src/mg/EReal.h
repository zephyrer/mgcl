/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGEReal_HH_
#define _MGEReal_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGEReal.h
//

class MGIfstream;
class MGOfstream;

/// MGEReal is extended real number, i.e., it includes minus infinite and 
/// plus infinite except ordinary real number.
class MGCLASS MGEReal {

public:

MGDECL friend MGEReal operator+ (double, const MGEReal&);
MGDECL friend MGEReal operator- (double, const MGEReal&);
MGDECL friend MGEReal operator* (double, const MGEReal&);
MGDECL friend MGEReal operator/ (double, const MGEReal&);
MGDECL friend bool operator== (double, const MGEReal&);
MGDECL friend bool operator!= (double, const MGEReal&);
MGDECL friend bool operator> (double, const MGEReal&);
MGDECL friend bool operator< (double, const MGEReal&);
MGDECL friend bool operator>= (double, const MGEReal&);
MGDECL friend bool operator<= (double, const MGEReal&);

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGEReal& );

////////////Constructor////////////

///Default Constructor
MGEReal(double val=0.0):m_value(val){;};

///infinite=-1 means minus_infinite,
///         +1 means plus_infinite.
MGEReal(MGINFINITE_TYPE infinite);

//Destructor
//	~MGEReal();	

////////////Operator overload.////////////

MGEReal operator+ (double) const;
MGEReal operator+ (const MGEReal&) const;

MGEReal& operator+= (double);
MGEReal& operator+= (const MGEReal&);

///Unary minus.
MGEReal operator- () const;

MGEReal operator- (double) const;
MGEReal operator- (const MGEReal&) const;

MGEReal& operator-= (double);
MGEReal& operator-= (const MGEReal&);

MGEReal operator* (double) const;
MGEReal operator* (const MGEReal&) const;

MGEReal& operator*= (double);
MGEReal& operator*= (const MGEReal&);

MGEReal operator/ (double) const;
MGEReal operator/ (const MGEReal&) const;

MGEReal& operator/= (double);
MGEReal& operator/= (const MGEReal&);

bool operator== (double t) const;
bool operator== (const MGEReal&) const;

bool operator!= (double t) const{return !((*this)==t);};
bool operator!= (const MGEReal& er2) const{return !((*this)==er2);};

bool operator> (double t) const;
bool operator> (const MGEReal&) const;

bool operator< (double t) const;
bool operator< (const MGEReal& er2) const{return er2>(*this);};

bool operator>= (const MGEReal& er2)const;
bool operator<= (const MGEReal& er2)const{return er2>=(*this);};
bool operator>= (double t) const;
bool operator<= (double t) const;

////////////Member Function////////////

///return -1 if minus_infinite(), 1 if plus_infinite(), else 0.
int infinite_coef()const;

bool equal_base(double t, double base)const;
bool equal_base(const MGEReal& t,double base)const;
bool finite()const{return infinite_coef()==0;};
bool infinite()const{return infinite_coef()!=0;};
void invert(){m_value*=-1.;};
bool minus_infinite()const{return (m_value<=(-mgInfiniteVal));};
bool plus_infinite()const{return (mgInfiniteVal<=m_value);};
void set_real(double val){m_value=val;};
void set_plus_infinite(){m_value=mgInfiniteVal+1.;};
void set_minus_infinite(){m_value=-mgInfiniteVal-1.;};
void set_zero(){m_value=0.;};
double value() const {return m_value;};

///Dump Functions.
///Calculate dump size
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );
	
////////////Member Data////////////

private:
	///int m_infinite;	-1: minus_infinite, 0:ordinary real, 1:plus_infinite
	double	m_value;///<when m_value<=-mgInfiniteVal, this is minus infinite.
					///<when m_value>=mgInfiniteVal, this is plus infinite.

friend class MGInterval;

};

/** @} */ // end of BASE group
#endif
