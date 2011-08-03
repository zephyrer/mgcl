/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLTexPoints_HH_
#define _mgTLTexPoints_HH_

#include "mg/MGCL.h"
#include <vector>
#include "mg/Position.h"

#if defined(MGCL_DLL)
#pragma warning(push)
#pragma warning(disable : 4231)
#pragma warning(disable : 4660)
MGTEMPLATE template class MGCLASS std::vector<MGPosition>;
#pragma warning(pop)
#endif

///mgTLTexPoints holds the vector of texture coordinates of mgTLPoints.
///the surface parameter (u,v).
class MGCLASS mgTLTexPoints{

public:

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLTexPoints& texPoints);

//////////// constructor ///////////////
mgTLTexPoints(){;};
mgTLTexPoints(
	size_t n			///<number of points.
);

////////// member function /////////////

MGPosition& operator[](size_t pos){return m_tex_positions[pos];};
const MGPosition& operator[](size_t pos) const{return m_tex_positions[pos];};

///get reference number.
int ref(size_t i){return m_refs[i];}

void set_texcoord(size_t i, double s, double t);
void set_texcoord(size_t i, const MGPosition &texc);
void set_texcoord_as_not_modify(size_t i, double s, double t);
void set_texcoord_as_not_modify(size_t i, const MGPosition &texc);
size_t size()const{return m_tex_positions.size();};
bool was_set_as_not_modify(size_t i)const{
	return m_refs[i]==-1;
}

private:
	std::vector<MGPosition> m_tex_positions;///<vector of texture coordinates
	std::vector<int> m_refs;///<m_refs[i] is the reference counter of m_tex_positions[i].
			///<m_refs[i]<0 means not to modify the texture data since it is set
			///<from the partner edge's data.
	///<double m_image_width;	//image width in world coordinates.
	///<double m_image_height;	//image height in world coordinates.

};

#endif
