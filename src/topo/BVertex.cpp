/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Point.h"
#include "topo/BVertex.h"
#include "topo/PVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGBVertex Class.
//MGBVertex is 0 manifold dimension binder cell, is an point.
//MGBVertex is a binder cell of MGPVertex, and has manifold dimension 0.

/////// Constructor ///////

//Copy constructor.
//MGBVertex(const MGBVertex&v);

//Fundamental constructor.
//Construct a BVertex from geometry geo of manifold dimension 0
//(MGPoint*, may be null).
//The constructor takes the ownership of geo.
MGBVertex::MGBVertex(MGGeometry* geo):MGCellNB(geo){;}

//Construct from MGPosition data.
MGBVertex::MGBVertex(const MGPosition& V):MGCellNB(new MGPoint(V)){;}

MGBVertex::~MGBVertex(){
	free_from_parent();//Free this cell from parent complex.
						//This is a Binder cel.
}

/////// operator overload///////

MGBVertex& MGBVertex::operator=(const MGBVertex& gel2){
	if(this==&gel2)
		return *this;
	MGCellNB::operator=(gel2);
	return *this;
}
MGBVertex& MGBVertex::operator=(const MGGel& gel2){
	const MGBVertex* gel2_is_this=dynamic_cast<const MGBVertex*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGBVertex::operator<(const MGBVertex& gel2)const{
	const MGPoint* pnt1=point();
	if(!pnt1)
		return true;

	const MGPoint* pnt2=gel2.point();
	if(!pnt2)
		return false;

	return (*pnt1)<(*pnt2);
}
bool MGBVertex::operator<(const MGGel& gel2)const{
	const MGBVertex* gel2_is_this=dynamic_cast<const MGBVertex*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

///////Member Function///////

//Make a clone of this(this is a binder), and set binder and member
//partner relation between the new binder and  the cell c.
MGBVertex* MGBVertex::clone_binder(const MGCellBase& c) const{
	assert(c.manifold_dimension()==0);

	MGBVertex* bvertex=new MGBVertex(*this);
	c.set_binder(*bvertex);
	return bvertex;
}

//Make a binder cell of this cell.
//Returned is the binder pointer newed.
//The binder has no geometry, only has binder and partner member relationship.
MGCellNB* MGBVertex::make_binder() const{
	MGCellNB* cell=binder();
	if(cell) return cell;
	MGBVertex* bvertex=new MGBVertex();
	set_binder(*bvertex);
	return bvertex;
}

//Make this cell's binder cell's extent expression.
//Returned is a MGGeometry pointer generated by new.
//When this cell does not have star cell, null pointer will be returned.
//make_binder_extent() only makes the expression, and does nothing to
//the topology structure.
MGGeometry* MGBVertex::make_binder_extent() const{
	const MGCellNB* cell=star();
	if(!cell) return 0;
	const MGFace* face=dynamic_cast<const MGFace*>(cell);
	make_extent();//Make sure that this has an extent expression.
	const MGPoint* vertx=point(); if(!vertx) return 0;
	MGPoint* pnt=new MGPoint(face->eval(vertx->position()));
	return pnt;
}

//Make sure that this has an extent expression.
//When this did not have an extent, make the extent from the partner
//member's parameter expression and the star cell.
//This must be a binder cell that has partner members that are
//boundaries. When this is not the case or this had an extent already,
//it does nothing.
void MGBVertex::make_extent() const{
	if(extent()) return;
	if(!is_bcell()) return;

	const MGCellBase* partner=m_partners[0];
	const MGCellNB* cell=partner->star();
	if(!cell) return;

	MGPosition pos(1);
	const MGPVertex* pv=dynamic_cast<const MGPVertex*>(partner);
	if(pv) pos(0)=pv->t();
	else{
		const MGBVertex* bv=dynamic_cast<const MGBVertex*>(partner);
		if(bv) pos=bv->point()->position();
		else return;
	}
	MGPoint* pnt=new MGPoint(cell->extent()->evaluate(pos));
	MGBVertex* thisB=const_cast<MGBVertex*>(this);
	thisB->set_extent(pnt);
}

//Obtain the i-th member partner vertex.
const MGPVertex* MGBVertex::member_partner_vertex(size_t i)const{
	return static_cast<const MGPVertex*>(member_partner(i));
}

// Output virtual function.
std::ostream& MGBVertex::out(std::ostream& ostrm) const{
	ostrm<<"<<BV=";
	MGCellNB::out(ostrm);
	ostrm<<"=BV>>";
	return ostrm;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
//This will be never invoked.
MGPosition MGBVertex::pick_closest(const MGStraight& sl)const{
	return point()->pick_closest(sl);
};

//Return extent data(i.e. MGPoint*).
const MGPoint* MGBVertex::point() const{
	return dynamic_cast<const MGPoint*>(extent());
}

//Get the position data of this vertex.
//Returns MGPoint data if this has extent. Otherwise obtain from partner's
//star edge by evaluating the edge's data.
MGPosition MGBVertex::position() const{
	const MGPoint* pnt=point();
	if(pnt)
		return pnt->position();

	const MGCellBase* partner=m_partners[0];
	const MGEdge* edg=dynamic_cast<const MGEdge*>(partner->star());
	if(!edg)
		return MGPosition();

	const MGPVertex* pv=dynamic_cast<const MGPVertex*>(partner);
	if(!pv)
		return MGPosition();

	return edg->eval(pv->t());
}