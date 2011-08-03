/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "Tl/TLPoints.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

size_t mgTLPoints::add(double u, double v){
	size_t id=m_positions.size();
	m_positions.push_back(MGPosition(u,v));
	return id;
}

ostream& operator<< (ostream& out, const mgTLPoints& tlPoints){
	out<<"TLPoints="<<(&tlPoints)<<" num="<<tlPoints.m_positions.size()<<endl;
	mgTLPoints::const_iterator iter = tlPoints.begin();
	size_t cnt = 0;
	for(; iter != tlPoints.end(); iter++){
		out<<cnt++<<" "<<*iter<<", ";
		if((cnt%4)==0)out<<endl;
	}
	out<<endl;
	return out;
}
