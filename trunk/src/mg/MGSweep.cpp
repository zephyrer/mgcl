/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/Straight.h"
#include "mg/SBRep.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//The sweep surface is defined as:
//lbrep(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSBRep::MGSBRep(
	const MGLBRep& lbrep,		//Sweep crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist)			//distance to end edge.
:MGSurface()
{
	size_t lBdim = lbrep.bdim();
	MGBPointSeq bpnt1 = lbrep.line_bcoef();		//曲線のB表現
	MGSPointSeq spnt1(lBdim, 2, 3);	//面のB表現

	//面B表現を作成する
	for(size_t i=0; i < lBdim; i++){
		//スイープ始点
		MGVector tmpSPnt = bpnt1(i);
		tmpSPnt += (uvec * start_dist);
		spnt1.store_at(i, 0, tmpSPnt);

		//スイープ終点
		MGVector tmpEPnt = bpnt1(i);
		tmpEPnt += (uvec * end_dist);
		spnt1.store_at(i, 1, tmpEPnt);
	}

	size_t sOrder = 2;
	size_t sBdim = 2;
	m_uknot = lbrep.knot_vector();
	m_vknot = MGKnotVector(sOrder, sBdim, 0.0, fabs(end_dist - start_dist));
	m_surface_bcoef = spnt1;
}

//The sweep surface is defined as:
//st(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSBRep::MGSBRep(
	const MGStraight& st,		//Sweep crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist)			//distance to end edge.
:MGSurface(){
	MGSPointSeq spnt1(2, 2, 3);	//面のB表現

	//直線の始終点を与えられたベクトル、長さで移動し、面B係数に入力する
	MGVector tmpSPnt = st.start_point();
	MGVector tmpEPnt = st.start_point();
	tmpSPnt += (uvec * start_dist);
	tmpEPnt += (uvec * end_dist);
	spnt1.store_at(0, 0, tmpSPnt);
	spnt1.store_at(0, 1, tmpEPnt);

	tmpSPnt = st.end_point();
	tmpEPnt = st.end_point();
	tmpSPnt += (uvec * start_dist);
	tmpEPnt += (uvec * end_dist);
	spnt1.store_at(1, 0, tmpSPnt);
	spnt1.store_at(1, 1, tmpEPnt);

	size_t sOrder = 2;
	size_t sBdim = 2;

	m_uknot = MGKnotVector(sOrder, sBdim, 0.0, st.param_span());
	m_vknot = MGKnotVector(sOrder, sBdim, 0.0, fabs(end_dist - start_dist));
	m_surface_bcoef = spnt1;
}

//The sweep surface is defined as:
//rlbrep(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGRSBRep::MGRSBRep(
	const MGRLBRep& rlbrep,		//Sweep crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist)			//distance to end edge.
:MGSurface(){
	size_t lBdim = rlbrep.bdim();
	MGBPointSeq rbpnt1 = rlbrep.homogeneous().line_bcoef();
	MGSPointSeq spnt1(lBdim, 2, 4);	//面のB表現(weight含むからsdim+1)

	//面B表現を作成する
	MGBPointSeq staBpnt = rbpnt1, endBpnt = rbpnt1;
	staBpnt.homogeneous_transform(uvec * start_dist);
	endBpnt.homogeneous_transform(uvec * end_dist);

	for(size_t i=0; i < lBdim; i++){
		spnt1.store_at(i, 0, staBpnt(i));
		spnt1.store_at(i, 1, endBpnt(i));
	}

	size_t sOrder = 2;
	size_t sBdim = 2;
	MGKnotVector uknot = rlbrep.knot_vector();
	MGKnotVector vknot = MGKnotVector(sOrder, sBdim, 0.0, fabs(end_dist - start_dist));

	m_surface = MGSBRep(spnt1, uknot, vknot);
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGLBRep::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{return new MGSBRep(*this, uvec, start_dist, end_dist);}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This straight(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGStraight::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{return new MGSBRep(*this, uvec, start_dist, end_dist);}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGRLBRep::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{return new MGRSBRep(*this, uvec, start_dist, end_dist);}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGEllipse::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{return new MGRSBRep((MGRLBRep)*this, uvec, start_dist, end_dist);}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGTrimmedCurve::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	MGCurve *tempCrv = clone();
	MGSurface *rtnSrf = tempCrv->sweep(uvec, start_dist, end_dist);
	delete tempCrv;
	return rtnSrf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGSurfCurve::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	MGLBRep tempCrv(*this);
	MGSurface *rtnSrf = tempCrv.sweep(uvec, start_dist, end_dist);
	return rtnSrf;
}
