/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/SurfCurve.h"

#include "topo/Complex.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Boundary.h"
#include "topo/Loop.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//--------------------------------------------
// MGTopologyのシリアライズ関数群
// 

//identify_type implements.
long MGComplex::identify_type()const{return MGCOMPLEX_TID;}
long MGLoop::identify_type()const{return MGLOOP_TID;}
long MGShell::identify_type()const{return MGSHELL_TID;}
long MGPVertex::identify_type()const{return MGPVERTEX_TID;}
long MGBVertex::identify_type()const{return MGBVERTEX_TID;}
long MGEdge::identify_type()const{return MGEDGE_TID;}
long MGFace::identify_type()const{return MGFACE_TID;}

//--------------------------------------------
// メンバデータ書き込み関数
void MGTopology::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
}
// メンバデータ読み出し関数
void MGTopology::ReadMembers(MGIfstream& buf){
	MGObject::ReadMembers(buf);
}

//--------------------------------------------
// メンバデータ書き込み関数
void MGCellBase::WriteMembers(MGOfstream& buf)const{
	MGTopology::WriteMembers(buf);
	buf.WritePointer(m_binder);
}
// メンバデータを読み出す関数
void MGCellBase::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGTopology::ReadMembers(buf);
	m_binder = static_cast<MGCellNB*>(buf.ReadPointer());
}

//--------------------------------------------
// メンバデータ書き込み関数
void MGCellNB::WriteMembers(MGOfstream& buf)const{
	MGCellBase::WriteMembers(buf);
	buf.WritePointer(m_parent_complex);
	buf.WritePointer(m_extent);
	size_t n=m_partners.size();
	buf<<n;
	for(size_t i=0; i<n; i++) buf.WritePointer(m_partners[i]);
}
// メンバデータを読み出す関数
void MGCellNB::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGCellBase::ReadMembers(buf);
	m_parent_complex=static_cast<MGComplex*>(buf.ReadPointer());
	m_extent=static_cast<MGGeometry*>(buf.ReadPointer());
	size_t n;
	buf>>n;//Number of partners.
	for(size_t i=0; i<n; i++)
		m_partners.push_back(static_cast<MGCellBase*>(buf.ReadPointer()));
}

//--------------------------------------------
// メンバデータ書き込み関数
void MGPVertex::WriteMembers(MGOfstream& buf)const{
	MGCellBase::WriteMembers(buf);
	buf.WritePointer(m_edge);
	buf<<m_t;
}
// メンバデータ読み出し関数
void MGPVertex::ReadMembers(MGIfstream& buf){
	MGCellBase::ReadMembers(buf);
	m_edge=static_cast<MGEdge*>(buf.ReadPointer());
	buf>>m_t;
}

//--------------------------------------------
// MGCellのシリアライズ関数群

// メンバデータを書き込む関数
void MGCell::WriteMembers(MGOfstream& buf)const{
	// 親クラスのデータ書き込み関数を呼んでおく。
	MGCellNB::WriteMembers(buf);
	m_box.dump(buf);
	size_t n=m_boundaries.size();
	buf<<n;
	for(size_t i=0; i<n; i++) buf.WritePointer(m_boundaries[i]);
	buf<<m_perror;
}

// メンバデータを読み出す関数
void MGCell::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGCellNB::ReadMembers(buf);
	m_box.restore(buf);
	size_t n;
	buf>>n;
	for(size_t i=0; i<n; i++)
		m_boundaries.push_back(static_cast<MGBoundary*>(buf.ReadPointer()));
	buf>> m_perror;
}

//--------------------------------------------
// MGBVertexのシリアライズ関数群
// 書き出しサイズを調べる関数

//Write Object's Member Data
void MGBVertex::WriteMembers(MGOfstream& buf) const{
	//親クラスのメンバデータの書き込み。
	MGCellNB::WriteMembers(buf);
}
//Read Object's member data.
void MGBVertex::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGCellNB::ReadMembers(buf);
}

//--------------------------------------------
//Write Object's Member Data
void MGEdge::WriteMembers(MGOfstream& buf) const{
	//親クラスのメンバデータの書き込み。
	//MGCellNB::WriteMembers(buf);
	MGCellBase::WriteMembers(buf);
	buf.WritePointer(m_parent_complex);
	const MGSurfCurve* scrv=dynamic_cast<const MGSurfCurve*>(m_extent);
	if(scrv) buf.WritePointer(0);
	else buf.WritePointer(m_extent);
	size_t n=m_partners.size();
	buf<<n;
	for(size_t i=0; i<n; i++) buf.WritePointer(m_partners[i]);
	//End of MGCellNB::WriteMembers(buf);
	
	buf.WritePointer(m_vertex[0]);
	buf.WritePointer(m_vertex[1]);
	m_box.dump(buf);
	buf<<m_perror;
}

//Read Object's member data.
void MGEdge::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGCellNB::ReadMembers(buf);
	m_vertex[0]=static_cast<MGPVertex*>(buf.ReadPointer());
	m_vertex[1]=static_cast<MGPVertex*>(buf.ReadPointer());
	m_box.restore(buf);
	buf>>m_perror;
}

//--------------------------------------------
// MGFaceのシリアライズ関数群

// メンバデータを書き込む関数
void MGFace::WriteMembers(MGOfstream& buf)const{
	//親クラスのメンバデータの書き込み。
	MGCell::WriteMembers(buf);
	m_box_param.dump(buf);
}
// メンバデータを読み出す関数
void MGFace::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGCell::ReadMembers(buf);
	m_box_param.restore(buf);
}

//--------------------------------------------
// MGComplexのシリアライズ関数群
// 

// メンバデータ書き込み関数
void MGComplex::WriteMembers(MGOfstream& buf)const{
	// 親クラスのデータ書き込み関数を呼んでおく。
	MGTopology::WriteMembers(buf);
	size_t n=m_pcells.size();
	buf<<n;
	MGComplex::const_pcellItr i=m_pcells.begin(), ie=m_pcells.end();
	for(; i!=ie; i++) buf.WritePointer(*i);
	n=m_bcells.size();
	buf<<n;
	MGComplex::const_bcellItr j=m_bcells.begin(), je=m_bcells.end();
	for(; j!=je; j++) buf.WritePointer(*j);
	m_box.dump(buf);
}

// メンバデータ読み出し関数
void MGComplex::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGTopology::ReadMembers(buf);

	size_t i,n;
	buf>>n;
	for(i=0; i<n; i++)
		m_pcells.push_back(static_cast<MGCellNB*>(buf.ReadPointer()));
	buf>>n;
	for(i=0; i<n; i++)
		m_bcells.push_back(static_cast<MGCellNB*>(buf.ReadPointer()));
	m_box.restore(buf);
}

//--------------------------------------------

// メンバデータを書き込む関数
void MGBoundary::WriteMembers(MGOfstream& buf)const{
	//親クラスのメンバデータの書き込み。
	MGComplex::WriteMembers(buf);
	buf.WritePointer(m_parent_cell);
}
// メンバデータを読み出す関数
void MGBoundary::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGComplex::ReadMembers(buf);
	m_parent_cell =static_cast<MGCell*>(buf.ReadPointer());
}

//--------------------------------------------
// MGLoopのシリアライズ関数群

// メンバデータを書き込む関数
void MGLoop::WriteMembers(MGOfstream& buf)const{
	//親クラスのメンバデータの書き込み。
	MGBoundary::WriteMembers(buf);
	buf << m_area << int(m_kind)
		<< m_perim_num_s << m_perim_num_e;
}
// メンバデータを読み出す関数
void MGLoop::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGBoundary::ReadMembers(buf);
	int kind;
	buf >> m_area >> kind
		>> m_perim_num_s >> m_perim_num_e;
	m_kind=LoopKind(kind);
}

//--------------------------------------------
// メンバデータを書き込む関数
void MGShell::WriteMembers(MGOfstream& buf)const{
	//親クラスのメンバデータの書き込み。
	MGBoundary::WriteMembers(buf);
}
// メンバデータを読み出す関数
void MGShell::ReadMembers(MGIfstream& buf){
	//親クラスのメンバデータの読み出し。
	MGBoundary::ReadMembers(buf);
}
