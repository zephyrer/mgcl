/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/LEPoint.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGShell Class.

///////Constructor////////

//Construct a shell of one face.
//Copy version of face.
MGShell::MGShell(const MGFace& face){
	MGFace* nf=new MGFace(face);
	nf->make_outer_boundary();
	append_pcell(nf);
}

//Construct a shell of one face.
//face is a pointer to a newed face, and the ownership
//will be transfered to MGShell.
MGShell::MGShell(MGFace* face){
	face->make_outer_boundary();
	append_pcell(face);
}

//Construct a shell of one face.
//Copy version of face.
MGShell::MGShell(const MGFSurface& face){
	MGFace* f=face.clone_as_face();
	f->make_outer_boundary();
	append_pcell(f);
}

//Construct a shell of one face.
//face is a pointer to newed face, and the ownership
//will be transfered to MGShell.
MGShell::MGShell(MGFSurface* face){
	MGFace* f=face->make_face();
	f->make_outer_boundary();
	append_pcell(f);
}

//Copy constructor with mapping.
//Binder cells of the pcells in loop will be registered in cmap
MGShell::MGShell(
	const MGShell& shell,	//original shell.
	MGCellMap& cmap)		//cellmap to register binder association.
//Binder cells of the pcells in loop will be registered in cmap.
:MGBoundary(shell, cmap){;}

//Fundamental constructor.
//Construct from boundary complex(i.e. MGLoop).
//This constructor takes the ownership of MGCell* in boundary.
MGShell::MGShell(
	std::list<MGCellNB*> boundaries
		//Boundary data of the super class MGBoundary(List of faces).
):MGBoundary(boundaries){;}

///////Destructor///////

///////operator overload///////
bool MGShell::operator<(const MGShell& gel2)const{
	return MGComplex::operator<(gel2);
}

bool MGShell::operator<(const MGGel& gel2)const{
	const MGShell* gel2_is_this=dynamic_cast<const MGShell*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGShell& MGShell::operator=(const MGShell& gel2){
	if(this==&gel2)
		return *this;

	MGBoundary::operator=(gel2);
	return *this;
}
MGShell& MGShell::operator=(const MGGel& gel2){
	const MGShell* gel2_is_this=dynamic_cast<const MGShell*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

///////Member Function///////

//Make a clone.
MGShell* MGShell::clone(MGCell& parent) const{
	MGShell* shl=clone();
	shl->set_parent(parent);
	return shl;
}
MGShell* MGShell::clone() const{
	MGShell* shl=new MGShell(*this);
	return shl;
}

//Make a clone.
//The forms that have cmap as an argumetnt is to register binder association.
//Binders of pcells in this shell will be registered in cmap.
//Returned is pointer of newed object, must be deleted.
//When parent is specified, clone's parent is set to the parent.
MGShell* MGShell::clone(MGCell& parent,MGCellMap& cmap) const{
	MGShell* shl=clone(cmap);
	shl->set_parent(parent);
	return shl;
}
MGShell* MGShell::clone(MGCellMap& cmap) const{
	MGShell* shl=new MGShell(*this,cmap);
	return shl;
}

//Make a clone that has not binders.
MGShell* MGShell::clone_without_binders(MGCell& parent) const{
	MGShell* shl=clone_without_binders();
	shl->set_parent(parent);
	return shl;
}
MGShell* MGShell::clone_without_binders() const{
	MGShell* shl=new MGShell();
	shl->copy_boundary_without_binders(*this);
	return shl;
}

//Exlude non free edges or opposite direction's edges from their vectors.
void exclude(std::vector<MGLEPoint>& pr1,
			 std::vector<MGLEPoint>& pr2){
	if(pr1.empty()) return;

	std::vector<MGLEPoint>::iterator i1=pr1.begin(), i1e=pr1.end(), i2e=pr2.end();
	//Exclude non free edges, or opposite direction's edges.
	do{
		i1e--, i2e--;
		if( (i1e->eval_star(1)%i2e->eval_star(1))>0. || !(i1e->edge()->is_free())
			|| !(i2e->edge()->is_free()) ){
			pr1.erase(i1e); pr2.erase(i2e);
		}
	}while(i1e!=i1);
}

//Get face pointer from its iterator in MGComplex of MGBoudarynD.
MGFace* MGShell::face(pcellItr i){
	return static_cast<MGFace*>(*i);
}
const MGFace* MGShell::face(const_pcellItr i)const{
	return static_cast<const MGFace*>(*i);
}

//Get the face pointer from its pcell id in MGComplex of MGBoudarynD.
MGFace* MGShell::face(size_t i){
	pcellItr j=pcell_begin();
	std::advance(j,i);
	return face(j);
}
const MGFace* MGShell::face(size_t i)const{
	const_pcellItr j=pcell_begin();
	std::advance(j,i);
	return face(j);
}

//Ensure this shell has binder vertices of all the model edges.
//The binder vertices are stored in bcells of this shell(in MGComplex).
void MGShell::ensure_BVertices_of_ModelEdges()const{
	int nfaces=number_of_faces();
	for(int i=0; i<nfaces; i++){
		const MGFace& fi=*(face(i));
		int nloops=fi.number_of_loops();
		for(int j=0; j<nloops; j++){
			const MGLoop& lj=*(fi.loop(j));
			if(lj.is_inactive())
				continue;
			int nedges=lj.number_of_edges();
			for(int k=0; k<nedges; k++){
				const MGEdge& ek=*(lj.edge(k));
				MGEdge* bek=ek.binder_edge();
				size_t id1, id2;
				if(ek.equal_direction_to_binder())
					id1=0;//to start of bek.
				else
					id1=1;//to the end of bek.
				const MGEdge& ek_pre=*(ek.pre_edge());
				MGEdge* bek_to_connect=ek_pre.binder_edge();
				if(ek_pre.equal_direction_to_binder())
					id2=1;
				else
					id2=0;
				bek->connect_at_id(id1,bek_to_connect,id2);
			}
		}
	}
}

//This is a proprietry routine for merge_at_common_edge, will merge
//loops at free common edges.
bool merge_loop(bool can_negate, MGLoop& lp1, MGLoop& lp2){
	std::vector<MGLEPoint> pr1,pr2;
	std::vector<double> br1,br2;
	if(!lp1.common(lp2,pr1,br1,pr2,br2)) return false;

	//for(int iii=0; iii<pr1.size();iii++){std::cout<<pr1[iii]<<","<<pr2[iii]<<std::endl;}
	if(can_negate && pr1[0].eval_star(1)%pr2[0].eval_star(1)>0.){
		//negate the face.
		//std::cout<<*(lp2.face());
		lp2.face()->negate();
		//std::cout<<*(lp2.face());
		lp1.common(lp2,pr1,br1,pr2,br2);
	}
	exclude(pr1,pr2);
	//Exclude pairs of edges that are not the same direction as the first edge.
	if(pr1.empty()) return false;
	//for(int iii=0; iii<pr1.size();iii++){std::cout<<pr1[iii]<<","<<pr2[iii]<<std::endl;}

	std::vector<MGEdge*> e1s=lp1.subdivide(pr1);
	std::vector<MGEdge*> e2s=lp2.subdivide(pr2);
	size_t ne=e1s.size();
	for(size_t i=0; i<ne; i++)
		e1s[i]->connect(*(e2s[i]));
	return true;
}

//Merge a face at free edges of this shell.
//Function's return value is
//   false: merge was not done because no common edges were found.
//   true: merge of the shell and the face was done.
//      Case of true includes that merge was done only partialy.
//      This case occurs when after some common edges are merged,
//      some common edges are found to have contradictionary direction.
//      Merge is done only for the first common edges found to have
//      same direction.
//When the function's return value is false, the face will not be added to
//the shell. Merge negates the input face direction when this shell's direction
//is not the same as the input face's one along the common edge(s).
//The second form is a pointer version. The face must be newed object
//and will be destructed after merge when function's return value is 0
//(merge was processed).
bool MGShell::merge_at_common_edge(const MGFace& face){
	MGFace* nf=new MGFace(face);
	bool merged=merge_at_common_edge(nf);
	if(!merged) delete nf;
	return merged;
}
bool MGShell::merge_at_common_edge(const MGFSurface& face){
	MGFSurface* nf=face.clone_fsurface();
	bool merged=merge_at_common_edge(nf);
	if(!merged) delete nf;
	return merged;
}
bool MGShell::merge_at_common_edge(MGFSurface* face){
	MGFace* f=dynamic_cast<MGFace*>(face);
	if(f)
		return merge_at_common_edge(f);

	MGSurface* srf=dynamic_cast<MGSurface*>(face);
	MGFace* f2=new MGFace(srf);
	if(merge_at_common_edge(f2))
		return true;
	f2->free_extent();
	delete f2;
	return false;
}

//This form is a pointer version. The face must be newed object
//and will be destructed after merge when function's return value is 0
//(merge was processed).
//This shell may be dummy shell. In this case, the face is added to the 1st
//face in the shell.
bool MGShell::merge_at_common_edge(MGFace* face){
	//std::cout<<(*this)<<std::endl;////*********
	size_t i1,i2, i1save,i2save; //id of inner boudary of this shell's loop and face.
	size_t nin1,nin2; //number of inner boudaries of this shell's loop and face.
	size_t j1,j2;	//inner boudary loop variable.

	face->make_outer_boundary();
	//std::cout<<(*face)<<std::endl; std::cout<<(face->outer_boundary())<<std::endl;/////
	MGLoop* out2=face->loop(size_t(0));
	nin2=face->number_of_inner_boundaries(i2save);

	pcellItr ci=pcell_begin(), ce=pcell_end();
	if(ci==ce){
		append_pcell(face);
		return true;
	}
	bool merged=false, merged2;
	for(;ci!=ce; ci++){
		MGFace* fi=dynamic_cast<MGFace*>(*ci);
		if(!fi) continue;
		if(fi==face) continue;

		//std::cout<<(*fi)<<std::endl;////*********
		// 1. outer1 versus outer2.
		MGLoop* out1=fi->loop(size_t(0));//std::cout<<(*out1);/////////
		merged2=merge_loop(!merged, *out1, *out2);
		if(merged2){
			merged=true;
			//std::cout<<(*out1);/////////
		}

		// 2. outer1 versus inner2.
		i2=i2save;
		for(j2=0; j2<nin2; j2++, i2++){
			merged2=merge_loop(!merged,*out1,*(face->loop(i2)));
			if(merged2) merged=true;
		}

		nin1=fi->number_of_inner_boundaries(i1save);
		i1=i1save;
		for(j1=0; j1<nin1; j1++, i1++){
		// 3. inner1 versus outer2.
			MGLoop* inner1=fi->loop(i1);
			merged2=merge_loop(!merged,*inner1, *out2);
			if(merged2) merged=true;

		// 4. inner1 versus inner2.
			i2=i2save;
			for(j2=0; j2<nin2; j2++, i2++){
				merged2=merge_loop(!merged,*inner1,*(face->loop(i2)));
				if(merged2) merged=true;
			}
		}
	}

	return merged;
}

//Debug Function
std::ostream& MGShell::out(std::ostream& ostrm) const{
	ostrm<<"<<MGShell="<<this<<"=";
	MGComplex::out(ostrm);
	ostrm<<"=MGShell>>"<<std::endl;
	return ostrm;
}

//Obtain boundary and main parameter lines of the FSurface.
//skeleton includes boundary() and inner parameter lines.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
MGPvector<MGCurve> MGShell::skeleton(int density) const{
	const_pcellItr ci=pcell_begin(), ce=pcell_end();
	MGPvector<MGCurve> crvs;
	for(;ci!=ce; ci++){
		MGFace* fi=dynamic_cast<MGFace*>(*ci);
		MGPvector<MGCurve> si=fi->skeleton(density);
		crvs.push_back(si);
	}
	return crvs;
}

//Obtain all the parameter curves at knots of u and v knot vector.
MGPvector<MGCurve> MGShell::skeleton_at_knots()const{
	const_pcellItr ci=pcell_begin(), ce=pcell_end();
	MGPvector<MGCurve> crvs;
	for(;ci!=ce; ci++){
		MGFace* fi=dynamic_cast<MGFace*>(*ci);
		MGPvector<MGCurve> si=fi->skeleton_at_knots();
		crvs.push_back(si);
	}
	return crvs;
}
