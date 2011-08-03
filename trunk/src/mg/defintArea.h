/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDefintArea_HH_
#define _MGDefintArea_HH_

/// @cond

#define mgdefintlen 608

class MGDefintArea{
private:
	MGDefintArea();
	MGDefintArea(const MGDefintArea&);
	MGDefintArea& operator=(const MGDefintArea&);

public:
	static MGDefintArea& instance();
	void init();

	double m_am[mgdefintlen*2]	/* was [608][2] */,
		m_a0[2], m_ap[mgdefintlen*2]	/* was [608][2] */,
		m_b0, m_bb[mgdefintlen];

	int m_nend;
	int m_npow;
	double m_eps0;
};

/// @endcond

#endif
