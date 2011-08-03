/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDrawFucn_HH_
#define _MGDrawFucn_HH_
#include "mg/MGCL.h"
#include "mgGL/Color.h"

/** @addtogroup ALGORITHM
 *  @{
 */

// MGDrawFunc.h
//

///This class is used for drawWire functions of all the MGObjects.
///MGDrawFunc includes four functions pointer. These functions are typically
///glBegin(xxxxx), glVertex3d(double,double,double), glEnd(),
///and a function to draw a point(DrawPoint below).
///Though glBegin needs an argument of the vertex mode, MGDrawFunc's first function
///does not. So, define the following type function. MGDrawFunc's first function is
///this BeginFunc:
///void BeginFunc(void){ glBegin(GL_LINE_STRIP);}
///DrawPoint is defined as:
///	void DrawPoint(double x,double y,double z);
class MGCLASS MGDrawFunc{

public:

typedef void (*BeginLineStripF)(void);
typedef void (*VertexF)(double, double, double);
typedef void (*EndF)(void);
typedef void (*DrawPoint)(double x,double y,double z);

////////////Constructor////////////

////////////Member Function////////////////

static BeginLineStripF bfunc(){ return m_bfunc;};
static VertexF vfunc(){ return m_vfunc;};
static EndF efunc(){ return m_efunc;};
static DrawPoint pfunc(){return m_pfunc;};

///Boundary color of the point.
const static MGColor& PColor1(){return m_PColor1;};

///Inner square color of the point.
const static MGColor& PColor2(){return m_PColor2;};

///Point rectangle width and height of below.
const static float point_size(){return m_point_size;};

static void set_bfunc(BeginLineStripF bfunc){ m_bfunc=bfunc;}
static void set_vfunc(VertexF vfunc){ m_vfunc=vfunc;}
static void set_efunc(EndF efunc){ m_efunc=efunc;}
static void set_pfunc(DrawPoint pfunc){ m_pfunc=pfunc;}

static void set_functions(
	BeginLineStripF bfunc,
	VertexF vfunc,
	EndF efunc,
	DrawPoint pfunc){
	m_bfunc=bfunc;
	m_vfunc=vfunc;
	m_efunc=efunc;
	m_pfunc=pfunc;
}

private:

////////////Member Data//////////
	static BeginLineStripF	m_bfunc;///<typically glBegin(GL_LINE_STRIP);
	static VertexF m_vfunc;		///<typically glVertex3d(double,double,double);
	static EndF	m_efunc;		///<typically glEnd();
	static DrawPoint m_pfunc;	///<typically DrawPoint(double x,double y,double z);
	const static MGColor& m_PColor1;///<Boundary color of the point.
	const static MGColor& m_PColor2;///<Inner square color of the point.
	const static float m_point_size;///<Point rectangle width and height of below.

};

/** @} */ // end of ALGORITHM group

#endif
