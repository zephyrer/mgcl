/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLFans_HH_
#define _mgTLFans_HH_

#include "mg/Pvector.h"
#include "Tl/TLFan.h"
#include "Tl/TLPoints.h"
#include "Tl/TLTriangle.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGPvector<mgTLFan>;
#pragma warning( pop )
#endif

///////////// mgTLFans /////////////

class MGPosition;
class mgTLEdges;
class mgTLTriangle;

///  private class for tessellation.
class MGCLASS mgTLFans{
public:
	typedef MGPvector<mgTLFan>::iterator iterator;
	typedef MGPvector<mgTLFan>::const_iterator const_iterator;
/////////////public member data//////////////////
	MGPvector<mgTLFan> m_fans;

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLFans& fans);

/////////////constructor////////////////////////
	mgTLFans(
		double uerror, double verror,///<Error along u and v in the parameter space.
		const mgTLPoints& tlPoints,	///<mgTLPoints of the tessellation input.
		mgTLTriangle& polygon);		///<このfansを作成するための元の多角形

/////////////operator overload////////////////////////
const mgTLFan* operator[](size_t i)const{return m_fans[i];};
mgTLFan* operator[](size_t i){return m_fans[i];};

/////////////member function////////////////////////

	const_iterator begin()const{return m_fans.begin();};
	const_iterator end()const{return m_fans.end();};
	iterator begin(){return m_fans.begin();};
	iterator end(){return m_fans.end();};

	///3点目の頂点(id of m_fans)を求める(頂点は使用点となる)
	size_t find3rdV(
		size_t		alpha,	///<エッジの始点(id of m_fans)
		size_t		beta,	///<エッジの終点(id of m_fans)
		size_t&		status	///<ステータス
	);

	///Fan生成に必要な変数の準備を行う
	///edgeStackにedgeのstackを積む
	void init(
		mgTLEdges& edgeStack
	);

	///Test if the edge(alpha, beta) is boundary or not.
	bool is_boundary(size_t alpha, size_t beta) const;

	///線分(v1,v2)が既存のedgeと交点があるかどうかを調べる
	///Function's return: true if had isect.
	bool has_isect(
		size_t 	v1,	///<線分1の点
		size_t 	v2	///<線分1の点
	)const;

	size_t pointID(size_t fanid)const{return m_polygon[fanid];};
	void push_back(mgTLFan* fan){m_fans.push_back(fan);};
	void push_back(mgTLFans& fans){m_fans.push_back(fans.m_fans);};

	///目的：中心点(center)がgammaのmgTLFanに対し､頂点addIndexを基準standIndexの
	///前（is_front=trueの時）､または後ろ（is_front=falseの時）追加する
	/*void push1V(
		size_t	gamma,		///中心点のインデックス
		size_t	addIndex,	///追加する頂点のインデックス
		size_t	standIndex,	///基準となる頂点のインデックス
		bool	is_front=false	///基準(standIndex)の前に追加するかどうかのフラグ
	);*/

	///目的：中心点(center)がalphaのmgTLFanに対し､頂点gammaを基準betaの
	///後ろに追加する
	void push1Vaft(
		size_t	alpha,	///<中心点のインデックス
		size_t	beta,	///<追加する頂点のインデックス
		size_t	gamma	///<基準となる頂点のインデックス
	);

	///目的：中心点(center)がalphaのmgTLFanに対し､頂点betaを基準gammaの
	///前に追加する
	void push1Vbefore(
		size_t	alpha,	///<中心点のインデックス
		size_t	beta,	///<追加する頂点のインデックス
		size_t	gamma	///<基準となる頂点のインデックス
	);

	///目的：中心点(center)がgammaで頂点(alpha,beta)のものを新規に作成する
	void push2V(
		size_t	gamma,	///<γのインデックス
		size_t	alpha,	///<αのインデックス
		size_t	beta	///<βのインデックス
	);

	///Set edge(alpha,j) as used.
	void set_edge_used(size_t alpha, size_t beta);

	///Get the number of fans included.
	size_t size() const{return m_polygon.size();};

	///check if vertex(alpha) is used or not.
	bool used(size_t alpha) const;

	///check if edge(alpha, beta) is used or not.
	bool used(size_t alpha, size_t beta) const;

	const MGPosition& uv(size_t i)const;

////////////////member data///////////////
private:
	double m_uerror;///<error along u in the parameter space of face or surface.
	double m_verror;///<error along v in the parameter space of face or surface.
	const mgTLPoints& m_points;	///<頂点列
	mgTLTriangle& m_polygon;///<このfansを作成するための元の多角形

};

#endif
