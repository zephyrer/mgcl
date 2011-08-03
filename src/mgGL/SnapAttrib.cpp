/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// SnapAttrib.cpp: MGSnapAttrib クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
#include "MGCLStdAfx.h"
#include "mgGL/SnapAttrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

/*	enum{
		BIT_END=0,
		BIT_KNOT=1,
		BIT_NEAR=2,
		BIT_VERTEX=3,
		BIT_CENTER=4,
		BIT_GRID=5,
		BIT_ORTHO=6
	};
*/

MGSnapAttrib::MGSnapAttrib()
:m_bitset(),m_dSnapApertureX(16.),m_dSnapApertureY(16.){;}

MGSnapAttrib::MGSnapAttrib(
	float apx, float apy,
	bool bEnd, bool bKnot, bool bNear, bool bVertex, bool bCenter, bool bGrid, bool bOrtho
):m_bitset(),m_dSnapApertureX(apx),m_dSnapApertureY(apy){
	m_bitset[BIT_END] = bEnd;
	m_bitset[BIT_KNOT] = bKnot;
	m_bitset[BIT_NEAR] = bNear;
	m_bitset[BIT_VERTEX] = bVertex;
	m_bitset[BIT_CENTER] = bCenter;
	m_bitset[BIT_GRID] = bGrid;
	m_bitset[BIT_ORTHO] = bOrtho;
}
MGSnapAttrib::MGSnapAttrib(
	float apx, float apy,
	const std::bitset<32>& bits
):m_bitset(bits),m_dSnapApertureX(apx),m_dSnapApertureY(apy){
}

//MGSnapAttrib::~MGSnapAttrib(){;}

//Text output to stream.
std::ostream& operator<< (std::ostream& ostrm, const MGSnapAttrib& atr){
	ostrm<<"SnapAttrib::"
		<< "End=" << atr.m_bitset[0] ;
	ostrm<< ", Knot=" << atr.m_bitset[1] << ", Near=" << atr.m_bitset[2];
	ostrm<< ", Vertex=" << atr.m_bitset[3] << ", Center=" << atr.m_bitset[4];
	ostrm<< ", Grid=" << atr.m_bitset[5] << ", Ortho=" << atr.m_bitset[6];
	ostrm<<", Aperture=("<< atr.m_dSnapApertureX<<","<< atr.m_dSnapApertureY<<")";
	return ostrm;
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGSnapAttrib& atr) {
	buf << atr.m_bitset.to_ulong();
	buf << atr.m_dSnapApertureX;
	buf << atr.m_dSnapApertureY;
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGSnapAttrib& atr) {
	unsigned long lbit;
	buf >> lbit;
	atr.m_bitset = std::bitset<32>(lbit);
	buf >> atr.m_dSnapApertureX;
	buf >> atr.m_dSnapApertureY;
	return buf;
}

/*
void MGSnapAttrib::ReadMembers(MGIfstream& buf){
	MGAttrib::ReadMembers(buf);

	unsigned long lbit;
	buf >> lbit;
	m_bitset = std::bitset<32>(lbit);
	buf >> m_dSnapApertureX;
	buf >> m_dSnapApertureY;

}

void MGSnapAttrib::WriteMembers(MGOfstream& buf) const{
	MGAttrib::WriteMembers(buf);

	buf << m_bitset.to_ulong();
	buf << m_dSnapApertureX;
	buf << m_dSnapApertureY;
}
*/
