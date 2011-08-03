/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgImageRect_HH_
#define _mgImageRect_HH_

#include <iosfwd>
#include <deque>
#include "mg/Position.h"
#include "mgGL/IRisect.h"

class MGBox;

/** @addtogroup DisplayHandling
 *  @{
 */

///mgImageRect is a proprietry class for texturemapping.
///mgImageRect holds all the necessary information for the triangulation
///of a retangle face parameter space.
class mgImageRect{

public:

typedef std::deque<mgIRisect>::iterator iterator;
typedef std::deque<mgIRisect>::const_iterator const_iterator;

friend std::ostream& operator<< (std::ostream& out, const mgImageRect& rect);

////////////Constructor/////////////
mgImageRect(
	const MGBox& box,///<The texture image box.
	const MGPosition st_tri[3],
	const MGPosition world_tri[3]
);

//////////// Member Function ///////////////

///Add vertices of the triangle to stpoly from the isect i to ip1.
///When ip1==m_isects.end(), ip1 means isects.begin().
void add_trim_point(
	iterator i,		///<m_isects' iterator.
	iterator ip1,	///<m_isects' iterator.
	std::vector<MGPosition>& stpoly
);

///Convert (s,t) coordinates to world coordinates using m_triangle and m_world_tri.
void convert_to_world(
	const MGPosition& st,	///<(s,t) coordinates
	MGPosition& world		///<world coordinates
)const;

///Convert (s,t) coordinates to world coordinates using m_triangle and m_world_tri.
void convert_to_world(
	const std::vector<MGPosition>& stpoly,
	std::vector<MGPosition>& worldpoly
)const;

///矩形交点からトリムポリゴンを作成し polygons に追加する
void createTrimPolygon(
	std::vector<MGPosition>& stpoly
);

///Test if (u,v) is in the triangle.
bool in_triangle(double u, double v)const;

///Get the start point parameter value (u,v) of the perimeter peri.
///Start point means in the order of the rectangle' anti-clockwise order.
MGPosition start_point(int peri)const;

///Return this rectangle's corner point i's (u,v).
MGPosition isect_uv(const mgIRisect& is)const;

private:
	double m_t[4];///Parameter range of the rectangle, that is,
			///m_t[]={vmin,umax,vmax,umin};

	const MGPosition* m_triangle[3];
	const MGPosition* m_world_tri[3];

	std::deque<mgIRisect> m_isects;	
	///<m_isects are intersections of the rectangle with the trinangle.
	///<In m_isects, intersections are sorted in the anti-clock-wise along
	///<the rectangle perimeter.
	///<On return from m_isects constructor,
	///<m_isects[2*j] is the point where the rectangle's perimeter is going
	///<into the trinangle and m_isects[2*j+1] is going out. These two points always makes a pair.

};

/** @} */ // end of DisplayHandling group
#endif
