/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFOuterCurve_HH_
#define _MGFOuterCurve_HH_

#include <iosfwd>
#include "mg/MGCL.h"
class MGLoop;

//
//Define MGFOuterCurve Class.

/** @addtogroup TOPORelated
 *  @{
 */

///MGFOuterCurve is to represent Face's outer boundary in combination of
///loops and perimeter curves. This is private class for Face and Loop.
class MGCLASS MGFOuterCurve{

public:

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGFOuterCurve& );

/////////Constructor/////////
MGFOuterCurve(){;};

///Constructor of perimeter curve.
MGFOuterCurve(
	size_t peri_id,			///<Perimeter id.
	double t1, double t2)	///<Parameter range of the perimeter curve.
	:m_peri_id(peri_id){m_t[0]=t1;m_t[1]=t2;};

///Constructor of loop.
MGFOuterCurve(
	const MGLoop* loop)
	:m_peri_id(-1), m_loop(loop){;};

/////////Operator oveload/////////

/////////Member function/////////

///Return true if this is perimeter curve.
bool is_perimeter_curve() const{return m_peri_id>=0;};

///Return true if this is loop.
bool is_loop() const{return m_peri_id<0;};

///Return loop pointer. This is valid only when is_loop is true.
const MGLoop* loop() const{return m_loop;};

///Return perimeter id. This is valid only when is_perimeter_curve is true.
size_t perimeter_id() const{return m_peri_id;};

///Return range of the perimeter curve. Valid only when is_perimeter.
void range(double& t0, double& t1)const{ t0=m_t[0]; t1=m_t[1];};

private:
	int m_peri_id;	///<perimeter number if >=0. When m_peri_id<0, 
					///<This is a loop and m_loop is active.
	const MGLoop* m_loop;///<When m_peri_id<0, m_loop is active and m_loop is 
					///<loop pointer of the Face.
	double m_t[2];	///<When m_peri_id>=0, m_t is active and indicates
					///<parameter range of the perimeter curve.

};

/** @} */ // end of TOPORelated group
#endif
