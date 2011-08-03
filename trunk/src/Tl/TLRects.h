/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLRects_HH_
#define _mgTLRects_HH_

#include "Tl/TLPoints.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<mgTLRect*>;
#pragma warning( pop )
#endif

class mgTLparameter;
class mgTLRect;
class mgTLpoints;
class mgTLTexPlane;

///mgTLRects is a proprietry class for Face tessellation.
///mgTLRects holds all the subdivided rectangles for the triangulation.
class MGCLASS mgTLRects{

public:

	typedef std::vector<mgTLRect*>::iterator RecIterator;
	typedef std::vector<mgTLRect*>::const_iterator const_RecIterator;
	typedef std::vector<mgTLRect*>::iterator iterator;
	typedef std::vector<mgTLRect*>::const_iterator const_iterator;
	std::vector<mgTLRect*> m_rects;	///<mgTLRect sequence ordered.

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLRects& rects);

//////////// constructor ///////////////
mgTLRects(
	mgTLparameter& param,	///<Tessellation parameter.
		///<Point((u,v) of the face) array will be returned into param.m_points.
		///<mgTLRect's m_Pid[.] is the id of this points.
			///<points[m_Pid[.]] is the MGPosition data (u,v) of the face.
	size_t divnum=4	///<Minimum number of the devided rectangles.
);

mgTLRects(const mgTLRects& rects2);	///Copy constructor.

////////// destructor ///////////////
~mgTLRects();

//////////// operator overload ////////////
mgTLRects& operator= (const mgTLRects& rects2);

////////// member function /////////////

///Add position data P into m_position.
///Function's return value i is the id of the stored position.
///m_position[i] will be P.
size_t add_point(double u, double v);

RecIterator begin(){return m_rects.begin();}
const_RecIterator begin()const{return m_rects.begin();}
const_RecIterator end()const{return m_rects.end();}
RecIterator end(){return m_rects.end();}

///Find mgTLRect that includes the surface parameter uv.
mgTLRect* find_rect(const MGPosition& uv);

///Test if surface of the rect is flat, that is, if within the tolerance.
///When m_param.texture() is true, test if rect is within the texture tolerance
///and if so, the texture plane will be pushed back to m_texture_planes.
bool is_flat(
	mgTLRect& rect,
	bool& direction///<  true: u-direction is more non flat.
					///< false: v-direction is more non flat.
);

mgTLparameter& parameter(){return *m_param;};

const mgTLPoints& points() const{return *m_points;}

void push_back(mgTLRect* rect){m_rects.push_back(rect);};
size_t size()const{return m_rects.size();};

const mgTLRect& rect(int i)const{return *(m_rects[i]);};

const std::vector<mgTLTexPlane*>& texture_planes()const{return m_texture_planes;};

private:
	mgTLparameter* m_param;	///<Only reference to mgTLparameter, not newed one for this TLRects.
	mgTLPoints*    m_points;///<Only reference to mgTLPoints, not newed one for this TLRects.
	std::vector<mgTLTexPlane*> m_texture_planes;///<vector of newed mgTLTexPlane to indicate
				///<texture mapping plane. This is valid only when m_param->texture() is true.
				///<Let i1=m_texture_planes[i]->index(), and i2=m_texture_planes[i+1]->index().
				///<Then rects of m_rects[i1] to m_rects[i2-1] are mapped onto the plane of
				///<m_texture_planes[i]->plane().

///init() builds the 1st rectangle information of the whole face.
///Function's return value is mgTLRect* of the 1st rectangle.
///init() also builds u=min and max, v=min and max mgTLisects of
///mgTLisectsList uList and vList in the m_param.
mgTLRect* init();

///Check if rect perimeter(that is rect's u=min,max or v=min,max parameter line)
///is within surface tolerance.
///If a rect  perimeter is not in surface tolerance,
///subdivide uvrect and push back to m_rects.
void rect_perimeter_subdivide(mgTLRect* uvrect);

///Check if surface perimeter is within curve tolerance.
///If the perimeter is not in curve tolerance, subdivide uvrect and push back to m_rects.
void surface_perimeter_subdivide(mgTLRect* uvrect);

///Check if ratio of u span and v span of uvrect is within maximum.
///If the ratio exceeds maximum, subdivide uvrect and push back to m_rects.
void ratio_subdivide(mgTLRect* uvrect);

friend class mgTLRect;

};

#endif
