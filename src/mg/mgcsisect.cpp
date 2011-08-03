/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "mg/Interval.h"
#include "mg/CSisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCSisect.cc
// MGCSisect 実装ファイル

//
// コンストラクタ
//
// 初期化なしで、交点を生成する。
MGCSisect::MGCSisect () {}

// 全てのコンポーネントを指定して交点を生成する
MGCSisect::MGCSisect(const MGPosition& isect, 
	double par, const MGPosition& uv, const MGCSRELATION rel)
	: m_ipoint(isect), m_t(par), m_uv(uv), m_rel(rel){;}

//////////Operator overload///////////

bool MGCSisect::operator< (const MGCSisect& csi)const{
	return m_t<csi.m_t;
}

bool MGCSisect::operator== (const MGCSisect& csi)const{
	return MGREqual(m_t,csi.m_t) && (m_uv==csi.m_uv);
}

//Ordering functions.
bool MGCSisect::operator< (const MGisect& is)const{
	const MGCSisect* cis=dynamic_cast<const MGCSisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}
bool MGCSisect::operator== (const MGisect& is)const{
	const MGCSisect* cis=dynamic_cast<const MGCSisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

/////// メンバ関数

//obtain the distance in parameter space.
void MGCSisect::distance(
	const MGCSisect& isect2,	//2nd isect.
	double& t, double& u, double& v
)const{
	t=fabs(m_t-isect2.m_t);
	u=fabs(m_uv[0]-isect2.m_uv[0]);
	v=fabs(m_uv[1]-isect2.m_uv[1]);
}

double MGCSisect::distance_square(const MGCSisect& isect2)
//Compute square of parameter space distance between this and
//isect2.
const{
	MGVector dif=point()-isect2.point();
	return dif%dif;
}

// Output virtual function.
std::ostream& MGCSisect::out(std::ostream& ostrm)const{
//	ostrm.setf ( ios::scientific, ios::floatfield );
//	ostrm.precision ( 10 );
	ostrm << "MGCSisect::"<<m_rel<<" m_ipoint="<<m_ipoint
		<<",m_t="<<m_t<<",m_uv="<<m_uv<<std::endl;
	return ostrm;
}
