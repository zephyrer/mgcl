/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGisect_HH_
#define _MGisect_HH_

/** @defgroup IsectContainer Intersection Containers
 *  Containers of intersection of MGObject subclasses, 
 *  MGCurve, MGFSurface, MGLoop, MGShell, etc.
 *  @{
 */

#include <iostream>
#include "mg/MGCL.h"
#include "mg/Position.h"
class MGObject;
class MGCurve;
class MGCCisect;
class MGCSisect;
class MGCFisect;
class MGSSisect;
class MGFFisect;

///Define MGisect Class.
///MGisect is an abstract class to provide unified interfaces to handle
///an intersection of two objects.
class MGCLASS MGisect{
friend class MGisects;

public:

////////Constructor/////////
	
////////Destructor/////////
virtual ~MGisect(){;};

////////Operator oveload////////

///Ordering functions.
virtual bool operator< (const MGisect& is)const=0;
virtual bool operator> (const MGisect& is)const{return is<(*this);};
virtual bool operator<= (const MGisect& is)const{return !(is<(*this));};
virtual bool operator>= (const MGisect& is)const{return !((*this)<is);};
virtual bool operator== (const MGisect& is)const=0;
virtual bool operator!= (const MGisect& is)const{return !operator==(is);};
virtual bool operator< (const MGCCisect& is)const=0;
virtual bool operator< (const MGCSisect& is)const=0;
virtual bool operator< (const MGCFisect& is)const=0;
virtual bool operator< (const MGSSisect& is)const=0;
virtual bool operator< (const MGFFisect& is)const=0;

////////Member function////////

///Exchange 1st and 2nd order of the parameter line representation.
virtual void exchange12()=0;

///Return the object of the intersection(world coordinates representation).
virtual const MGObject& isect()const=0;

///Return the 1st object's parameter value of the intersection.
///*****This function is valid only when manifold_dimension()==0.
virtual MGPosition isect0_param1()const{return MGPosition();};

///Return the 2nd object's parameter value of the intersection.
///*****This function is valid only when manifold_dimension()==0.
virtual MGPosition isect0_param2()const{return MGPosition();};

///Return the 1st object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
virtual const MGCurve* isect1_param1()const{return 0;};

///Return the 2nd object's parameter value curve of the intersection.
///*****This function is valid only when manifold_dimension()==1.
virtual const MGCurve* isect1_param2()const{return 0;};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
virtual int manifold_dimension()const=0;

/// Output virtual function.
virtual std::ostream& out(std::ostream& ostrm)const=0;

private:

	///Get the object1 pointer.
	virtual const MGObject* object1(const MGObject* obj)const{return obj;};

	///Get the object2 pointer.
	virtual const MGObject* object2(const MGObject* obj)const{return obj;};

};

///Debug Function
inline std::ostream& operator<< (std::ostream& ostrm, const MGisect& is){
	is.out(ostrm);
	return ostrm;
}

/** @} */ // end of IsectContainer group
#endif
