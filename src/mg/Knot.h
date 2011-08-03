/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGKnot_HH_
#define _MGKnot_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGKnot.h
//

class MGIfstream;
class MGOfstream;

/// Defines knot value and its multiplicity.
class MGCLASS MGKnot{

public:

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGKnot& );

////////////Constructor////////////
MGKnot():m_multiplicity(0){;};	///Default Constructor

///*****This is the fundamental constructor.*****
MGKnot(	double t,			///Knot value
		int multiplicity);	///Multiplicity

//Destructor
//	~MGKnot();	

////////////Operator overload.////////////

bool operator< (const MGKnot& kt2)const{return m_value<kt2.m_value;};
bool operator> (const MGKnot& kt2)const{return kt2<(*this);};
bool operator<= (const MGKnot& kt2)const{return !(kt2<(*this));};
bool operator>= (const MGKnot& kt2)const{return !((*this)<kt2);};
bool operator== (const MGKnot& kt2)const;
bool operator!= (const MGKnot& kt2)const{return !operator==(kt2);};

////////////Member Function////////////

double value() const {return m_value;}
int multiplicity() const {return m_multiplicity;}

///Dump Functions.
///Calculate dump size
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

////////////Member Data////////////

private:
	double	m_value;			///<Knot value.
	int		m_multiplicity;		///<multiplicity of the value.

};

/** @} */ // end of BASE group
#endif
