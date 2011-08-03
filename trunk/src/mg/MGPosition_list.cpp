/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/Plane.h"
#include "mg/SSisect.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

#define mg_PLerror 7.5e-4

// MGPosition_list defines singly linked list of MGPosition.

//	RWTValSlist<MGPosition> *m_list; // List of position values. 

// Constructor

/*#define mg_SSallowance 1000.
#define mg_CSallowance 500.
#define mg_CCallowance 200.
#define mg_Sallowance  400.
#define mg_Callowance  250.
//10.*10. (square of 10.) is error allowed to remove.
//If within this allowance, position data will be removed or will not stored.
*/

//Construct MGPosition_list by replacing each element P of list by Q,
//where, Q=MGPosition(P.sdim(),P,start1,start2).
MGPosition_list::MGPosition_list
	(const MGPosition_list& lst,	//Original MGPosition_list 
		size_t start1,	//Start position to retrieve of elements P of list.
		size_t start2)	//Start position to store of elements Q of
						//new MGPosition_list.
{
	const_iterator i;
	for(i=lst.begin(); i!=lst.end(); i++){
		MGPosition Q((*i).sdim(),(*i),start1,start2);
		push_back(Q);
	}
}

//Copy Constructor.

// Destructor.

// Operator overload.

//Assignment.

// Member Function.

//Add parameter pair uvuv of surface srf1 and sur2 to the end(append==true)
//or the beginning of list if uvuv is not included in the list already.
//Function's return value is true if appended, false if not.
bool MGPosition_list::add(
	bool append, 
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
){
// Adds the MGPosition to the end of the list.
	MGEReal span[4];
	if(!empty()){
		MGBox s1r=srf1.param_range(), s2r=srf2.param_range();
		span[0]=s1r[0].length(); span[1]=s1r[1].length();
		span[2]=s2r[0].length(); span[3]=s2r[1].length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
	const_iterator i;
	for(i=begin(); i!=end(); i++){
		double a;
		size_t j;
		for(j=0;j<4;j++){
			a=(*i).ref(j)-uvuv.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==4) {MGTolerance::set_rc_zero(rcSave);return false;}
	}
	if(append) push_back(uvuv); else push_front(uvuv);
	MGTolerance::set_rc_zero(rcSave);
	return true;
}

//Function's return value is true if appended, false if not.
bool MGPosition_list::append(const MGPosition& position){
// Appends the position to the end of the list.
	push_back(position);
	return true;
}

//Add parameter uv of surface to the end of list if uv is not
//included in the list already.
//Function's return value is true if appended, false if not.
bool MGPosition_list::append(
	const MGFSurface& srf,
	const MGPosition& uv){
	MGEReal span[2];
	if(!empty()){
		MGBox s1r=srf.param_range();
		span[0]=s1r[0].length(); span[1]=s1r[1].length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
// Adds the MGPosition to the end of the list.
	const_iterator i;
	for(i=begin(); i!=end(); i++){
		double a;
		size_t j;
		for(j=0;j<2;j++){
			a=(*i).ref(j)-uv.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==2) {MGTolerance::set_rc_zero(rcSave);return false;}
	}
	push_back(uv);
	MGTolerance::set_rc_zero(rcSave);
	return true;
}

//Add parameter uv of surface to the end of list if uv is not
//included in the list already.
void MGPosition_list::append(
	const MGFSurface& srf,
	const MGPosition_list& lst){
	// Adds the MGPosition to the end of the list.
	const_iterator i;
	for(i=lst.begin(); i!=lst.end(); i++) append(srf,(*i));
}

//Add parameter pair p of crv1 and 2 to the end of list if p is not
//included in the list already.
//Function's return value is true if appended, false if not.
bool MGPosition_list::append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition& p)
{
	MGEReal span[2];
	if(!empty()){
		span[0]=crv1.param_range().length();
		span[1]=crv2.param_range().length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
// Adds the MGPosition to the end of the list.
	const_iterator i;
	for(i=begin(); i!=end(); i++){
		double a;
		size_t j;
		for(j=0;j<2;j++){
			a=(*i).ref(j)-p.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==2) {MGTolerance::set_rc_zero(rcSave);return false;}
	}
	push_back(p);
	MGTolerance::set_rc_zero(rcSave);
	return true;
}

//Add parameter pair p of crv1 and 2 to the end of list if p is not
//included in the list already.
void MGPosition_list::append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition_list& lst
){
	// Adds the MGPosition to the end of the list.
	const_iterator i;
	for(i=lst.begin(); i!=lst.end(); i++) append(crv1,crv2,(*i));
}

//Add parameter pair p of crv and surface to the end of list if p is not
//included in the list already.
//Function's return value is true if appended, false if not.
bool MGPosition_list::append(
	const MGCurve& crv,
	const MGFSurface& srf,
	const MGPosition& p)
{
	MGEReal span[3];
	if(!empty()){
		span[0]=crv.param_range().length();
		MGBox s1r=srf.param_range();
		span[1]=s1r[0].length(); span[2]=s1r[1].length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
// Adds the MGPosition to the end of the list.
	const_iterator i;
	for(i=begin(); i!=end(); i++){
		double a;
		size_t j;
		for(j=0;j<3;j++){
			a=(*i).ref(j)-p.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==3) {MGTolerance::set_rc_zero(rcSave);return false;}
	}
	push_back(p);
	MGTolerance::set_rc_zero(rcSave);
	return true;
}

//Add parameter pair uvuv of surface srf1 and sur2 to the end of list
void MGPosition_list::append(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition_list& lst
){
	// Adds the MGPosition to the end of the list.
	const_iterator i=lst.begin(), ie=lst.end();
	for(; i!=ie; i++)
		append(srf1,srf2,(*i));
}

//Test if one of the points in the list is included in
//the box. If so, return id of the point. The id is of the first point
//found in the list.
//Function's return value is true if a point is in the box, and false
//if no points are in the box.
bool MGPosition_list::in(const MGBox& box, iterator& id,size_t n){
	iterator i, ids=begin(), ide=end();
	if(ids==ide)
		return false;

	if(!n)
		n=box.sdim();
	MGBox boxt(n,box);
	MGPosition midP=box.mid();
	double min_len=-1.;
	i=ide;
	do{
		i--;
		MGPosition Pi(n,*i);
		if(boxt>>Pi){
			MGVector diff=Pi-midP;
			double len2=diff%diff;
			if(min_len<0. || min_len>len2){
				min_len=len2;
				id=i;
			}
		}
	}while(i!=ids);
	return min_len>0.;
}

//Test if one of the points in the list is included in
//the box. If so, return id of the point. The id is of the first point
//found in the list.
//Function's return value is true if a point is in the box, and false
//if no points are in the box.
bool MGPosition_list::in(const MGBox& box, const_iterator& id,size_t n)const{
	const_iterator i, ids=begin(), ide=end();
	if(ids==ide)
		return false;

	if(!n)
		n=box.sdim();
	MGBox boxt(n,box);
	MGPosition midP=box.mid();
	double min_len=-1.;
	i=ide;
	do{
		i--;
		MGPosition Pi(n,*i);
		if(boxt>>Pi){
			MGVector diff=Pi-midP;
			double len2=diff%diff;
			if(min_len<0. || min_len>len2){
				min_len=len2;
				id=i;
			}
		}
	}while(i!=ids);
	return min_len>0.;
}

//Remove position in the list that is the same point as P.
//When distace of the two point is within error, they are regarded as same.
//Function's return value is true if removal was done, and false if
//no same point was found and removal was not performed.
int MGPosition_list::remove(
	double error,			//square of error allowed to regard same.
	const MGPosition& P,	//Position to remove.
	size_t n)				//Number of space dimension to check.
							//If n==0, all the coordinates are checked.
{
	size_t num_removed=0;
// Remove the same point as uvuv from the list.
	iterator i=begin(), inext;
	double err;
	while(i!=end()){
		inext=i;inext++;
		if(n) { MGVector dif(n,(*i) - P); err=dif%dif;}
		else { MGVector dif=(*i) - P; err=dif%dif;}
		if(err<=error){ removeAt(i); ++num_removed;}
		i=inext;
	}
	return num_removed;
}

//Remove parameter pair uv of surface srf1 from the list
// if uv is included in the list.
int MGPosition_list::remove(
	const MGFSurface& srf,
	const MGPosition& uv)
{
	MGEReal span[2];
	MGBox s1r=srf.param_range();
	span[0]=s1r[0].length(); span[1]=s1r[1].length();
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
// Adds the MGPosition to the end of the list.
	iterator i=begin(), inext;
	size_t num_removed=0;
	while(i!=end()){
		inext=i;inext++;
		double a;
		size_t j;
		for(j=0;j<2;j++){
			a=(*i).ref(j)-uv.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==2) { removeAt(i); ++num_removed;}
		i=inext;
	}
	MGTolerance::set_rc_zero(rcSave);
	return num_removed;
}

//Remove parameter pair uvuv of surface srf1 and sur2 from the list
// if uvuv is included in the list.
int MGPosition_list::remove(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv){
	MGEReal span[4];
	if(!empty()){
		MGBox s1r=srf1.param_range(), s2r=srf2.param_range();
		span[0]=s1r[0].length(); span[1]=s1r[1].length();
		span[2]=s2r[0].length(); span[3]=s2r[1].length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
	iterator i=begin(), inext;
	size_t num_removed=0;
	while(i!=end()){
		inext=i;inext++;
		double a;
		size_t j;
		for(j=0;j<4;j++){
			a=(*i).ref(j)-uvuv.ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==4) { removeAt(i); ++num_removed;}
		i=inext;
	}
	MGTolerance::set_rc_zero(rcSave);
	return num_removed;
}

//Remove elements in this list(parameter pair uvuv of surface srf1 and sur2)
//from the list if elements are within the tolerance.
//Space dimension of the elements(uvuv) is 4. And uvuv(0) and uvuv(1) are srf1's
// parameter (u,v), uvuv(2) and uvuv(3) are srf2's parameter (u,v).
//Function's return value is number of removed elements.
int MGPosition_list::remove_uvuv(
	const MGFSurface& srf1,
	const MGFSurface& srf2){
	size_t n=size();
	if(n<=2) return 0;

	MGEReal span[4];
	if(!empty()){
		MGBox s1r=srf1.param_range(), s2r=srf2.param_range();
		span[0]=s1r[0].length(); span[1]=s1r[1].length();
		span[2]=s2r[0].length(); span[3]=s2r[1].length();
	}
	double rcSave=MGTolerance::set_rc_zero(mg_PLerror);
	iterator i=begin(), inext;
	inext=i; inext++;
	size_t num_removed=0;
	bool aft=true;//Flag to judge if aft elements should be removed.
				//true:aft element, false:before element.
	for(size_t count=2; count<=n; count++){
		size_t j;
		for(j=0;j<4;j++){
			double a=(*i).ref(j)-(*inext).ref(j);
			if(!MGRZero2(a, span[j])) break;
		}
		if(j==4){
			if(aft){
				iterator isave=inext; isave++;
				removeAt(inext);
				inext=isave;
			}else{
				removeAt(i);
				i=inext;
				inext++;
			}
			++num_removed;
		}else{
			//‚±‚±‚ÅŒðü‚Ì“r’†‚Ì“_‚©‚Ç‚¤‚©‚Ì”»’è‚ð‚·‚éB
			aft=!aft;
		}
	}
	MGTolerance::set_rc_zero(rcSave);
	return num_removed;
}

MGPosition MGPosition_list::removeAt(iterator i){
//Remove the position and return the position. If i is no valid, 
// behavior is undefined.
	MGPosition p=*i;
	erase(i);
	return p;
}

MGPosition MGPosition_list::removeFirst(){
//Remove the first position int the list and return the position.
//If i is not valid, behavior is undefined.
	MGPosition p=front();
	pop_front();
	return p;
}

MGPosition MGPosition_list::removeLast(){
//Remove the last position in the list and return the position.
//If i is not valid, behavior is undefined.
	MGPosition p=back();
	pop_back();
	return p;
}

//Remove uvuv that is on the ssi.
//Function's return value is the number of points deleted.
size_t MGPosition_list::removeOn(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGSSisect& ssi){
	iterator i=begin(), inext;
	size_t num_removed=0;
	const MGCurve& worldline=ssi.line();
	const MGCurve& uvline1=ssi.param1();
	const MGCurve& uvline2=ssi.param2();
	while(i!=end()){
		inext=i;inext++;
		double u1=(*i)[0], v1=(*i)[1];
		double u2=(*i)[2], v2=(*i)[3];
		MGPosition P((srf1.eval(u1,v1)+srf2.eval(u2,v2))*.5);
		double t;
		if(worldline.on(P,t)){
			double error=MGTolerance::set_wc_zero(srf1.param_error());
			if(uvline1.on(MGPosition(u1,v1),t)){
				MGTolerance::set_wc_zero(srf2.param_error());
				if(uvline2.on(MGPosition(u2,v2),t)){
					removeAt(i); num_removed++;
				}
			}
			MGTolerance::set_wc_zero(error);
		}
		i=inext;
	}
	return num_removed;
}

//reverse the ordering of the elements in the list.
void MGPosition_list::reverse_order(){
	size_t n=size();
	size_t n2=n/2;
	for(size_t i=0; i<n2; i++){
		MGPosition a=removeFirst(), b=removeLast();
		push_front(b); push_back(a);
	}
}

class MGPosition_list_compare:public std::greater<MGPosition>{
public:
	MGPosition_list_compare(size_t id, double center[2], const MGPosition& P1);
	bool operator()(const MGPosition* pos1, const MGPosition* pos2)const;
	double m_center[2];
	double m_V[2];
	size_t m_id;
};
MGPosition_list_compare::MGPosition_list_compare(
	size_t id, double center[2], const MGPosition& P1
	):m_id(id){
	m_center[0]=center[0]; m_center[1]=center[1];
	m_V[0]=P1[0]-center[0]; m_V[1]=P1[1]-center[1]; 
}
bool MGPosition_list_compare::operator()(
	const MGPosition* pos1,
	const MGPosition* pos2)const{
	size_t idp1=m_id+1;
	double a10=pos1->ref(m_id)-m_center[0], a11=pos1->ref(idp1)-m_center[1];
	double a20=pos2->ref(m_id)-m_center[0], a21=pos2->ref(idp1)-m_center[1];
	double len1=a10*m_V[0]+a11*m_V[1];
	double len2=a20*m_V[0]+a21*m_V[1];
	return len1>len2;
}

//Sort the positions in the list according to the surface parameter space
//ordering. Positions (uvuv(id),uvuv(id+1)) in the list is a parameter
//position of a face.
void MGPosition_list::sort_uv_space(size_t id){
	const size_t n=size();
	if(n<=2) return;
	size_t idp1=id+1, j;
	double center[2]={0.,0.};
	std::vector<MGPosition*> posv(n);
	iterator i, is=begin(), iend=end();
	for(i=is,j=0;j<n; j++, i++){
		posv[j]=&(*i);
		center[0]+=(*i)[id];
		center[1]+=(*i)[idp1];
	}
	center[0]/=double(n); center[1]/=double(n);
	MGPosition_list_compare comp(id, center,*is);
//	cout<<"Before sort::"<<(*this);
	std::sort(posv.begin(), posv.end(),comp);
	MGPosition_list list2;
	for(j=0; j<n; j++) list2.append(*(posv[j]));
	(*this)=list2;
//	cout<<"After sort::"<<(*this)<<endl;
}
