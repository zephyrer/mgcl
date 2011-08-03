/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/defintArea.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//***************************************************** 
// Generate points and weights of the DE formula for DEFINT
// in double precision.
// ***************************************************** *
void MGDefintArea::init(){
	if(m_nend>0) return;
	m_nend=mgdefintlen;
	int nendm1=m_nend-1;
	m_npow=6;
	m_eps0=1e-32;
    const double a9 = .9999999999999998;

//     --- start computation of points and weights --- 
    double ph = atan(1.) * 2.;
    m_a0[0] = 0.;
    m_a0[1] = 1.;
    m_b0 = ph;

    double eh = exp(.0078125);
    double en = 1.;

    for(int i = 1; i <= m_nend; ++i){
		en = eh * en;
		double eni = 1. / en;
		double sh = (en - eni) * .5;
		double ch = (en + eni) * .5;
		double exs = exp(ph * sh);
		double exsi = 1. / exs;
		double chsi = 2. / (exs + exsi);
		m_ap[i - 1] = (exs - exsi) * .5 * chsi;
		if (m_ap[i - 1] >= a9) m_ap[i - 1] = a9;
		m_ap[i + nendm1] = exsi * chsi;
		m_am[i - 1] = -m_ap[i - 1];
		m_am[i + nendm1] = -m_ap[i + nendm1];
		// Computing 2nd power
		m_bb[i - 1] = ph * ch * (chsi * chsi);
    }

	m_ap[m_nend+nendm1] = m_ap[nendm1*2];
	m_am[m_nend+nendm1] = m_am[nendm1*2];
}

MGDefintArea::MGDefintArea() :
	m_eps0(1e-32),
	m_nend(-1),
	m_npow(-1),
	m_b0(0.)
{
	m_a0[0] = m_a0[1] = 0.;
	std::fill_n(m_am, mgdefintlen*2, 0.);
	std::fill_n(m_ap, mgdefintlen*2, 0.);
	std::fill_n(m_bb, mgdefintlen,   0.);
}

MGDefintArea& MGDefintArea::instance(){
	static MGDefintArea theInst;
	return theInst;
}