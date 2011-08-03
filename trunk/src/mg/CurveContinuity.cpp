/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// CurveContinuity.cpp
//
// geometric continuity module

#include "MGCLStdAfx.h"
#include "mg/CurveContinuity.h"
#include "mg/Curve.h"
#include "mg/Tolerance.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

// Measures geometric continuity for two curves.
// returns true if calculation is successful.
MGCurveContinuity::MGCurveContinuity(
	const MGCurve& curve1,
	const MGCurve& curve2
){
	// 1. decide endpoint.
	double t1s=curve1.param_s(), t1e=curve1.param_e();
	MGPosition P1S=curve1.eval(t1s), P1E=curve1.eval(t1e);

	double t2s=curve2.param_s(), t2e=curve2.param_e();
	MGPosition P2S=curve2.eval(t2s), P2E=curve2.eval(t2e);

	double d1s2s=P1S.distance(P2S);
	double d1e2s=P1E.distance(P2S);
	double d1s2e=P1S.distance(P2E);
	double d1e2e=P1E.distance(P2E);

	if(d1s2s<=d1e2s){
		if(d1s2s<=d1s2e){
			if(d1s2s<=d1e2e){
				m_dist=d1s2s;
				m_param1=t1s; m_param2=t2s;
			}else{
				m_dist=d1e2e;
				m_param1=t1e; m_param2=t2e;
			}
		}else{
			if(d1s2e<=d1e2e){
				m_dist=d1s2e;
				m_param1=t1s; m_param2=t2e;
			}else{
				m_dist=d1e2e;
				m_param1=t1e; m_param2=t2e;
			}
		}
	}else{
		if(d1e2s<=d1s2e){
			if(d1e2s<=d1e2e){
				m_dist=d1e2s;
				m_param1=t1e; m_param2=t2s;
			}else{
				m_dist=d1e2e;
				m_param1=t1e; m_param2=t2e;
			}
		}else{
			if(d1s2e<=d1e2e){
				m_dist=d1s2e;
				m_param1=t1s; m_param2=t2e;
			}else{
				m_dist=d1e2e;
				m_param1=t1e; m_param2=t2e;
			}
		}
	}

	m_continuity = MGAZero(m_dist) ? G0 : DISCONT;

	m_P1=curve1.eval(m_param1);
	m_P2=curve2.eval(m_param2);
	MGUnit_vector b1, b2;
	double tor1, tor2;
	curve1.Frenet_frame(m_param1, m_tan1, m_normal1, b1, m_curvature1, tor1);
	curve2.Frenet_frame(m_param2, m_tan2, m_normal2, b2, m_curvature2, tor2);
	m_tandiff = m_tan1.angle(m_tan2);
	m_normaldiff=m_normal1.angle(m_normal2);

	if(m_continuity==G0){
		// 2. check whether G1 or not
		if(m_tandiff <= MGTolerance::angle_zero()){
			m_continuity = G1;

			// 3. finally, check whether G2 or not
			if(MGMZero(m_curvature1) || MGMZero(m_curvature2)){
				if(MGMZero(m_curvature1) && MGMZero(m_curvature2)){
					m_continuity = G2;
				}
			}else{
				double radius1=1./m_curvature1;
				double radius2=1./m_curvature2;
				if(MGREqual(radius1,radius2) && m_normaldiff<=MGTolerance::angle_zero()){
					m_continuity = G2;
				}
			}
		}
	}
}
