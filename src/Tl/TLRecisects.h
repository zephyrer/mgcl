/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLRecisects_HH_
#define _mgTLRecisects_HH_

#include <vector>
#include "Tl/TLRecisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<mgTLRecisect>;
#pragma warning( pop )
#endif

///mgTLRecisects is a proprietry class for Face tessellation.
///mgTLRecisects holds all the intersections of the rectangle with a loop.
///In mgTLRecisects, intersections are sorted in the anti-clock-wise along
///the rectangle perimter.
///On return from mgTLRects constructor,
///m_Recisects[2*j] is the point that is going into the rectangle and
///m_Recisects[2*j+1] is going out. These two points always makes a pair.
class MGCLASS mgTLRecisects{

public:

	typedef std::vector<mgTLRecisect>::iterator ReciItr;
	typedef std::vector<mgTLRecisect>::const_iterator CReciItr;
	std::vector<mgTLRecisect> m_Recisects;	
	///<m_Recisects are intersections of the rectangle with a loop.
	///<In m_Recisects, intersections are sorted in the anti-clock-wise along
	///<the rectangle perimter.
	///<On return from mgTLRects constructor,
	///<m_Recisects[2*j] is the point that is going into the rectangle and
	///<m_Recisects[2*j+1] is going out. These two points always makes a pair.

//////////// constructor ///////////////
mgTLRecisects(){;};
mgTLRecisects(size_t n):m_Recisects(n){;};

//////////// Operator overload ///////////////
mgTLRecisect& operator[](size_t i){ return m_Recisects[i];};
const mgTLRecisect& operator[](size_t i)const{ return m_Recisects[i];};

//////////// Member Function ///////////////

///Return the begin and end of m_Recisects.
ReciItr begin(){ return m_Recisects.begin();};
ReciItr end(){ return m_Recisects.end();};
CReciItr begin()const{ return m_Recisects.begin();};
CReciItr end()const{ return m_Recisects.end();};

///Erase the member at ReciItr.
void clear(){m_Recisects.clear();};

///Erase the member at ReciItr.
void erase(ReciItr i){m_Recisects.erase(i);};

///Push back a Recisect to the end.
void push_back(const mgTLRecisect& ri){m_Recisects.push_back(ri);};

///Resize the length of m_Recisects.
void resize(size_t n){m_Recisects.resize(n);};

///Get the size of m_Recisects.
size_t size()const{return m_Recisects.size();};

friend class mgTLRect;

};

#endif
