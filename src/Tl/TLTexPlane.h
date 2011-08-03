/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLTexPlane_HH_
#define _mgTLTexPlane_HH_

#include <iosfwd>

#include "mg/Plane.h"

///////////// mgTLTexPlane /////////////

//////////// private class for texture mapping ////////////

///mgTLTexPlane holds a plane to texture-map.
///The plane's origin is the center of the plane and the x-axis is the x-axis of
///the rectangle to hold the texture mapping.
///The width and the height is the rectangle data to cover texture mapping area.
class MGCLASS mgTLTexPlane{
private:

	MGPlane m_plane;	///<the plane fro teh texture mapping.
	double 	m_width;	///<The width of the texture mapping rectangle.
	double 	m_height;	///<The height of the texture mapping rectangle.
	int m_index;		///<the index of the mgTLRects to indicate the start of the
						///<mgTLRect for this mgTLTexPlane.
						///<From m_rects[i] are mgTLRect for this mgTLTexPlane.

public:

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLTexPlane& texplane);

///////////////Constructor///////////////

	mgTLTexPlane():m_width(0.),m_height(0.),m_index(-1){;};
	mgTLTexPlane(const MGPlane& plane,double width,double height,int index)
		:m_plane(plane),m_width(width),m_height(height),m_index(index){;};

/////////////member functions////////////

	const MGPlane& plane()const{return m_plane;};
	double width()const{return m_width;};
	double height()const{return m_height;};
	int index()const{return m_index;};
};

#endif
