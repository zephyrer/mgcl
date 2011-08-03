/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgIRisect_HH_
#define _mgIRisect_HH_

///@cond

///mgIRisect is a proprietry class for the class mgImageRect.
///mgIRisect represents one intesection of an mgImageRect perimeter and u=const(or v=const)
///line. When intersection with u=const, m_t is v value, or vise versa.
class mgIRisect{

public:
friend std::ostream& operator<< (std::ostream& out, const mgIRisect& isect);

mgIRisect():m_perimeter(0),m_t(0.),m_tri_id(0){;};
mgIRisect(int perimeter, double t, int tri_id)
:m_perimeter(perimeter),m_t(t),m_tri_id(tri_id){;};

/////////Operator oveload/////////

///Comparison with mgIRisect.
bool operator> (const mgIRisect& is)const;
bool operator< (const mgIRisect& is)const{return is>(*this);};
bool operator>= (const mgIRisect& is)const{return !(*this<is);};
bool operator<= (const mgIRisect& is)const{return !(*this>is);};
bool operator== (const mgIRisect& is)const;
bool operator!= (const mgIRisect& is)const{return !((*this)==is);};

int perimeter()const{return m_perimeter;};
double t()const{return m_t;};
int tri_id()const{return m_tri_id;};

private:
	int m_perimeter;///<Perimeter number of the rectangle, from 0 to 3.
	double m_t;	///<u or v value of the rectangle according to the perimeter number.
			///<When perimeter num=0, 2: u value and =1,3:v value.

	int m_tri_id;	///<Triangle's perimeter id.
};

///@endcond

#endif
