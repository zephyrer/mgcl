/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTopology_HH_
#define _MGTopology_HH_

#include "mg/Object.h"

class MGBox;
class MGCell;
class MGCellBase;
class MGComplex;
class MGOfstream;
class MGIfstream;
class MGVector;
class MGMatrix;
class MGTransf;

//Define MGTopology Class.

/** @defgroup TOPO Topology (sub) classes
 *  MGTopology is top abstract class for MGBVertex, MGPVertex, MGEdge, MGLoop,
 *  MGFace, MGShell.
 *  @{
 */

///MGTopology is an abstract class which represents a whole Topology,
///Complex, Cell, and Boundary.
///*****MGTopology must not have data members since this is an interface.*****
class MGCLASS MGTopology: public MGObject{

public:
///////// Constructor /////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGTopology();

///////// Destructor /////////
virtual ~MGTopology();

///Assignment.
///When the leaf object of this and obj2 are not equal, this assignment
///does nothing.
virtual MGTopology& operator=(const MGTopology& gel2);

/////////Member Function/////////

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
virtual MGisects intersection(const MGObject& obj2)const;
virtual MGisects intersection(const MGCurve& obj2)const;
virtual MGisects intersection(const MGFSurface& obj2)const;
virtual MGisects intersection(const MGSurface& obj2)const;
virtual MGisects intersection(const MGFace& obj2)const;
virtual MGisects intersection(const MGShell& obj2)const;

///Return MGTopology pointer if this MGGel is an MGTopology, else return null.
MGTopology* topology(){return this;};
const MGTopology* topology()const{return this;};

virtual std::string whoami()const{return "Topology";};

protected:

///Write Object into file stream
///It returns current streampos (=PID)
virtual void WriteMembers(MGOfstream& buf) const;

///Read Object from file stream
///It returns current streampos (=PID)
virtual void ReadMembers(MGIfstream& buf);

};

/** @} */ // end of TOPO group
#endif
