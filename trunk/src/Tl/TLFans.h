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
		mgTLTriangle& polygon);		///<����fans���쐬���邽�߂̌��̑��p�`

/////////////operator overload////////////////////////
const mgTLFan* operator[](size_t i)const{return m_fans[i];};
mgTLFan* operator[](size_t i){return m_fans[i];};

/////////////member function////////////////////////

	const_iterator begin()const{return m_fans.begin();};
	const_iterator end()const{return m_fans.end();};
	iterator begin(){return m_fans.begin();};
	iterator end(){return m_fans.end();};

	///3�_�ڂ̒��_(id of m_fans)�����߂�(���_�͎g�p�_�ƂȂ�)
	size_t find3rdV(
		size_t		alpha,	///<�G�b�W�̎n�_(id of m_fans)
		size_t		beta,	///<�G�b�W�̏I�_(id of m_fans)
		size_t&		status	///<�X�e�[�^�X
	);

	///Fan�����ɕK�v�ȕϐ��̏������s��
	///edgeStack��edge��stack��ς�
	void init(
		mgTLEdges& edgeStack
	);

	///Test if the edge(alpha, beta) is boundary or not.
	bool is_boundary(size_t alpha, size_t beta) const;

	///����(v1,v2)��������edge�ƌ�_�����邩�ǂ����𒲂ׂ�
	///Function's return: true if had isect.
	bool has_isect(
		size_t 	v1,	///<����1�̓_
		size_t 	v2	///<����1�̓_
	)const;

	size_t pointID(size_t fanid)const{return m_polygon[fanid];};
	void push_back(mgTLFan* fan){m_fans.push_back(fan);};
	void push_back(mgTLFans& fans){m_fans.push_back(fans.m_fans);};

	///�ړI�F���S�_(center)��gamma��mgTLFan�ɑ΂�����_addIndex���standIndex��
	///�O�iis_front=true�̎��j��܂��͌��iis_front=false�̎��j�ǉ�����
	/*void push1V(
		size_t	gamma,		///���S�_�̃C���f�b�N�X
		size_t	addIndex,	///�ǉ����钸�_�̃C���f�b�N�X
		size_t	standIndex,	///��ƂȂ钸�_�̃C���f�b�N�X
		bool	is_front=false	///�(standIndex)�̑O�ɒǉ����邩�ǂ����̃t���O
	);*/

	///�ړI�F���S�_(center)��alpha��mgTLFan�ɑ΂�����_gamma���beta��
	///���ɒǉ�����
	void push1Vaft(
		size_t	alpha,	///<���S�_�̃C���f�b�N�X
		size_t	beta,	///<�ǉ����钸�_�̃C���f�b�N�X
		size_t	gamma	///<��ƂȂ钸�_�̃C���f�b�N�X
	);

	///�ړI�F���S�_(center)��alpha��mgTLFan�ɑ΂�����_beta���gamma��
	///�O�ɒǉ�����
	void push1Vbefore(
		size_t	alpha,	///<���S�_�̃C���f�b�N�X
		size_t	beta,	///<�ǉ����钸�_�̃C���f�b�N�X
		size_t	gamma	///<��ƂȂ钸�_�̃C���f�b�N�X
	);

	///�ړI�F���S�_(center)��gamma�Œ��_(alpha,beta)�̂��̂�V�K�ɍ쐬����
	void push2V(
		size_t	gamma,	///<���̃C���f�b�N�X
		size_t	alpha,	///<���̃C���f�b�N�X
		size_t	beta	///<���̃C���f�b�N�X
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
	const mgTLPoints& m_points;	///<���_��
	mgTLTriangle& m_polygon;///<����fans���쐬���邽�߂̌��̑��p�`

};

#endif
