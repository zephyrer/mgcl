/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// SnapAttrib.h: MGSnapAttrib クラスのインターフェイス
#if !defined(AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_)
#define AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif /// _MSC_VER > 1000

#include "mg/MGCL.h"
#include <bitset>
#include <iosfwd>

class MGIfstream;
class MGOfstream;

///Defines Snap attributes.
///Snap means snap to designated location when inputting positional data
///from Cursor(generally mouse position).
///Currently MGCL supports the following enum's snap kind.
class MGCLASS MGSnapAttrib{
public:

	enum{
		BIT_END=0,
		BIT_KNOT=1,
		BIT_NEAR=2,
		BIT_VERTEX=3,
		BIT_CENTER=4,
		BIT_GRID=5,
		BIT_ORTHO=6
	};

///Text output to stream.
MGDECL friend std::ostream& operator<< (std::ostream& ostrm, const MGSnapAttrib& atr);

/// Serialization.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGSnapAttrib& atr);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGSnapAttrib& atr);

	MGSnapAttrib();
	MGSnapAttrib(float apx, float apy,
		bool bEnd=false, bool bKnot=false, bool bNear=false, bool bVertex=false,
		bool bCenter=false, bool bGrid=false, bool bOrtho=false
	);
	MGSnapAttrib(
		float apx, float apy,
		const std::bitset<32>& bits
	);

///	virtual ~MGSnapAttrib();

public:

	/// Attribute reference
	float getSnapApertureX()const{return m_dSnapApertureX;};
	void setSnapApertureX(float dApertureX){m_dSnapApertureX = dApertureX;};
	void setSnapApertureX(double dApertureX){m_dSnapApertureX = float(dApertureX);};
	float getSnapApertureY()const{return m_dSnapApertureY;};
	void setSnapApertureY(float dApertureY){m_dSnapApertureY = dApertureY;};
	void setSnapApertureY(double dApertureY){m_dSnapApertureY = float(dApertureY);};

	bool getEnd()  const {return m_bitset[BIT_END];};
	void setEnd(bool bEnd) {m_bitset[BIT_END] = bEnd;};

	bool getKnot() const {return m_bitset[BIT_KNOT];};
	void setKnot(bool bKnot) {m_bitset[BIT_KNOT] = bKnot;};

	bool getNear() const {return m_bitset[BIT_NEAR];};
	void setNear(bool bNear) {m_bitset[BIT_NEAR] = bNear;};

	bool getVertex() const {return m_bitset[BIT_VERTEX];};
	void setVertex(bool bVertex) {m_bitset[BIT_VERTEX] = bVertex;};

	bool getCenter() const { return m_bitset[BIT_CENTER];};
	void setCenter(bool bCenter) {m_bitset[BIT_CENTER] = bCenter;};

	bool getGrid() const { return m_bitset[BIT_GRID];};
	void setGrid(bool bGrid) {m_bitset[BIT_GRID] = bGrid;};

	bool getOrtho() const{ return m_bitset[BIT_ORTHO];}
	void setOrtho(bool bOrtho){ m_bitset[BIT_ORTHO] = bOrtho;}

	///////////// Proxy interfaces.//////////

	///Test if any of the attributes are set on. If on, return true.
	bool any() const{ return m_bitset.any();}

	///Test if none of the attributes is set on. If none is on, return true.
	bool none() const{ return m_bitset.none();}

	///Get the number of attributes that are on.
	size_t count() const{ return m_bitset.count();}

	///clear all the attrib to false;
	void clear_attrib(){m_bitset.reset();};

///member data
private:
	std::bitset<32> m_bitset;
	float m_dSnapApertureX;
	float m_dSnapApertureY;

};

#endif // !defined(AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_)
