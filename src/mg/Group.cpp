/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Group.h"
#include "mg/Attrib.h"
#include "mgGL/Appearance.h"
#include "mgGL/Context.h"
#include "mg/GelFactory.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGroup Class.
//MGGroup is a class which constains MGGel elements.
//MGGroup provides functions:
//(1) Serialization of MGGel elements.
//

//////////Constructor//////////

//Void constructor(初期化なしでオブジェクトを作成する。)
MGGroup::MGGroup(){;}
MGGroup::MGGroup(MGGel* gel){m_gels.push_back(gel);}

//Construct MGGroup from the file made by MGOfstream.
//error contains error flag as:
//=0: open succeeded.
//=1: file not found, or could not be opened.
//=2: file found, but, the format is not MGCL format.
//When error!=0, MGGroup contains no MGGel's.
MGGroup::MGGroup(const char* file, int& error){
	MGIfstream infile;
	error=infile.open(file,false);
	if(error) return;
	infile>>(*this);
}

///////////////operator overloaded//////////////

//Assignment.
//When the leaf objects of this and gel2 are not equal, this assignment
//does nothing.
MGGroup& MGGroup::operator=(const MGGroup& gel2){
	if(this==&gel2)
		return *this;

	MGAttribedGel::operator=(gel2);
	m_gels.clear();
	const_iterator i=gel2.begin(), ie=gel2.end();	
	for(; i!=ie; i++){
		m_gels.push_back((*i)->clone());
	}
	return *this;
}
MGGroup& MGGroup::operator=(const MGGel& gel2){
	const MGGroup* gel2_is_this=dynamic_cast<const MGGroup*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGGroup::operator<(const MGGroup& gel2)const{return size()<gel2.size();};
bool MGGroup::operator<(const MGGel& gel2)const{
	const MGGroup* gel2_is_this=dynamic_cast<const MGGroup*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//////////Member Function//////////

//Get the MGAppearance pointer in this group. If not defined, null will be
//returned.
MGAppearance* MGGroup::appearance(){
	if(!size())
		return 0;
	MGGel* gel=front();
	MGAppearance* app=dynamic_cast<MGAppearance*>(gel);
	if(app)
		return app;		
	MGContext* cnt=dynamic_cast<MGContext*>(gel);
	if(cnt)
		return cnt->appearance();		
	return 0;
}
const MGAppearance* MGGroup::appearance()const{
	MGGroup* grp=const_cast<MGGroup*>(this);
	return grp->appearance();
}

//Get the box of the group.
//If no objects were included in this group, null box will be returned.
MGBox MGGroup::box()const{
	MGBox bx;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(*i);
		if(obj){
			const MGBox& bxobj=obj->box();
			if(bxobj.finite())
				bx|=bxobj;
			continue;
		}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(*i);
		if(grp){
			bx|=grp->box();
			continue;
		}
	}
	return bx;
}

//Generate copied gel of this gel.
//Returned is a newed object. User must delete the object.
MGGroup* MGGroup::clone()const{
	MGGroup* gel=new MGGroup;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		gel->m_gels.push_back((**i).clone());
	}
	return gel;
}

//Get the MGContext pointer stored in this group. If not defined, null will be
//returned.
MGContext* MGGroup::context(){
	if(!size()) return 0;
	return dynamic_cast<MGContext*>(front());
}
const MGContext* MGGroup::context()const{
	if(!size()) return 0;
	return dynamic_cast<const MGContext*>(front());
}

//////display member function.
void MGGroup::display_arrows()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; ++i)
		(**i).display_arrows();
}
void MGGroup::display_break_points()const{
	MGGroup::const_iterator i=begin(), ie=end();
	for(; i!=ie; ++i)
		(**i).display_break_points();
}
void MGGroup::display_control_polygon()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; ++i)
		(**i).display_control_polygon();
}
void MGGroup::display_curvatures(
	double	scale,	//scaling of the graph.
	int		density,//densitiy of the graph.
	bool	use_radius//true:radius display, false:curvature display.
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; ++i)
		(**i).display_curvatures(scale,density,use_radius);
}

//make this group has appearance and get the MGAppearance pointer.
MGAppearance* MGGroup::ensure_appearance(){
	MGAppearance* app=appearance();
	if(!app){
		app=new MGAppearance();
		push_appearance(app);
	}
	return app;
}

//Find the position of the gel in the gel list and the group pointer
//which includes the gel. Searching will be done into the member group gel
//of this list.
MGGroup::const_iterator MGGroup::find(
	const MGGel* gel, const MGGroup*& grp
)const{
	std::vector<const MGGroup*> grps;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		if(*i==gel){ grp=this; return i;}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(*i);
		if(grp) grps.push_back(grp);
	}
	size_t ngrp=grps.size();
	for(size_t j=0; j<ngrp; j++){
		const MGGroup* gelj=grps[j];
		const_iterator pos=gelj->find(gel,grp);
		if(grp) return pos;
	}
	grp=0;
	return end();
}
MGGroup::iterator MGGroup::find(
	MGGel* gel, MGGroup*& grp
){
	std::vector<MGGroup*> grps;
	iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		if(*i==gel){ grp=this; return i;}
		MGGroup* grp=dynamic_cast<MGGroup*>(*i);
		if(grp) grps.push_back(grp);
	}
	size_t ngrp=grps.size();
	for(size_t j=0; j<ngrp; j++){
		MGGroup* gelj=grps[j];
		iterator pos=gelj->find(gel,grp);
		if(grp) return pos;
	}
	grp=0;
	return end();
}

//Test if this gel includes an object.
//Function's return value is the 1st object found in the gel list of this
//if this includes an object. Otherwise null will be returned.
const MGObject* MGGroup::includes_object()const{
	const_reverse_iterator i=rbegin(), ie=rend();	
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(*i);
		if(obj)
			return obj;
		const MGGroup* grp=dynamic_cast<const MGGroup*>(*i);
		if(grp){
			obj=grp->includes_object();
			if(obj)
				return obj;
		}
	}
	return 0;
}
MGObject* MGGroup::includes_object(){
	return
		const_cast<MGObject*>(const_cast<const MGGroup*>(this)->includes_object());
}

//Make an MGOfstream file which contains this group.
//The file generated by make_file can be retrieved by the constructor
//MGGroup(const char* file, int error);
//Function's return value is:
//=0: the file is successfully made.
//=1: file could not be opened.
int MGGroup::make_file(const char* file){
	MGOfstream outfile;
	int error=outfile.open(file,false);
	if(error)
		return error;
	outfile<<(*this);
	return 0;
}

//Make a display list of only call lists of this group.
//Return is the display list name.
//When gels in this group are included in gels_to_delete, the display list
//will not be made.
size_t MGGroup::make_only_call_list(
	bool no_color,	//if true, color attribute will be neglected.
	const std::vector<const MGGel*>* gels_to_delete
)const{
	size_t name=dlist_name();

	//make the wire mode display list to call list.
	int mode=0;
	if(no_color)
		mode=1;
	make_only_call_list_sub(name,mode,gels_to_delete);
	if(no_color)
		return name;

	//make the shading mode display list to call list.
	make_only_call_list_sub(name,2,gels_to_delete);
	return name;
}

// Output virtual function.
std::ostream& MGGroup::out(std::ostream& ostrm) const{
	ostrm<<"MGGroup::number of gels="<<size()<<std::endl;
	const_iterator i=begin(), ie=end();	
	for(size_t j=0; i!=ie; i++, j++){
		ostrm<<"gel"<<j<<":"<<(**i)<<std::endl;
	}
	return ostrm;
}

//Get the number of objects included in thie group.
int MGGroup::num_of_objects() const{
	int num=0;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(*i);
		if(obj){ num++; continue;}
		const MGGroup* grp=dynamic_cast<const MGGroup*>(*i);
		if(grp){num+=grp->num_of_objects(); continue;}
	}
	return num;
}

void MGGroup::push_appearance(MGAppearance* appr){
	if(!size()){
		m_gels.push_back(appr);
		return;
	}
	MGContext* ctx=dynamic_cast<MGContext*>(front());
	if(ctx){
		ctx->set_appearance(appr);
		return;
	}
	MGAppearance* app2=dynamic_cast<MGAppearance*>(front());
	if(app2) pop_front();
	m_gels.push_front(appr);
}
void MGGroup::push_context(MGContext* cntx){
	if(!size()){
		m_gels.push_back(cntx);
		return;
	}
	MGContext* ctx2=dynamic_cast<MGContext*>(front());
	if(ctx2){
		pop_front();
		m_gels.push_front(cntx);
		return;
	}
	MGAppearance* app2=dynamic_cast<MGAppearance*>(front());
	if(app2){
		release(begin());
		MGAppearance* app3=cntx->appearance();
		if(!app3) cntx->set_appearance(app2);
		else delete app2;
	}
	m_gels.push_front(cntx);
}

//push elements in objs at the end. All of the object pointers are
//transfered to this. On return, objs will have no object pointer in it.
void MGGroup::push_back(MGPvector<MGObject>& objs){
	size_t n=objs.size();
	for(size_t i=0; i<n; i++){
		MGObject* obj=objs.release(i);
		push_back(obj);
	}
	objs.clear();
}
void MGGroup::push_back(MGPlist<MGObject>& objs){
	size_t n=objs.size();
	for(size_t i=0; i<n; i++){
		MGObject* obj=objs.front();
		push_back(objs.front());
		objs.pop_front();
	}
}
void MGGroup::push_back(MGPvector<MGGel>& gels){
	size_t n=gels.size();
	for(size_t i=0; i<n; i++){
		MGGel* gel=gels.release(i);
		push_back(gel);
	}
	gels.clear();
}
void MGGroup::push_back(MGPlist<MGGel>& gels){
	size_t n=gels.size();
	for(size_t i=0; i<n; i++){
		MGGel* gel=gels.front();
		push_back(gel);
		gels.pop_front();
	}
}

// push element x at the end.
void MGGroup::push_back(MGGel* x){
	if(!size()){
		m_gels.push_back(x);
		return;
	}
	MGAppearance* app=dynamic_cast<MGAppearance*>(x);
	if(app){ push_appearance(app); return;}
	MGContext* ctx=dynamic_cast<MGContext*>(x);
	if(ctx){ push_context(ctx); return;}
	m_gels.push_back(x);
}

// push element x at the first.
void MGGroup::push_front(MGGel* x){
	if(!size()){
		m_gels.push_front(x);
		return;
	}
	MGAppearance* app=dynamic_cast<MGAppearance*>(x);
	if(app){ push_appearance(app); return;}
	MGContext* ctx=dynamic_cast<MGContext*>(x);
	if(ctx){ push_context(ctx); return;}

	iterator first=begin();
	app=dynamic_cast<MGAppearance*>(*first);
	if(app){
		first++; m_gels.insert(first,x);
		return;
	}
	ctx=dynamic_cast<MGContext*>(*first);
	if(ctx){
		first++; m_gels.insert(first,x);
		return;
	}
	m_gels.push_front(x);
};

//Remove the MGAppearance of this MGAttribedGel.
void MGGroup::remove_appearance(){
	if(!size()) return;

	MGGel* gel=front();
	MGAppearance* app=dynamic_cast<MGAppearance*>(gel);
	if(app){
		erase(begin());
		return;		
	}
	MGContext* cnt=dynamic_cast<MGContext*>(gel);
	if(cnt) cnt->set_appearance(0);		
}

//Read all member data.
void MGGroup::ReadMembers(MGIfstream& buf){
	MGGel::ReadMembers(buf);
	size_t n;
	buf>>n;
	for(size_t i=0; i<n; i++) push_back(buf.ReadPointer());
}
//Write all member data
void MGGroup::WriteMembers(MGOfstream& buf)const{
	MGGel::WriteMembers(buf);
	size_t n=size();
	buf<<n;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++) buf.WritePointer(*i);
}

//Transform the gel by the argument.
void MGGroup::transform(const MGVector& v)//translation
{
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++)
		(**i).transform(v);

}
void MGGroup::transform(double scale)//scaling.
{
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++)
		(**i).transform(scale);

}
void MGGroup::transform(const MGMatrix& mat)//matrix transformation.
{
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++)
		(**i).transform(mat);

}
void MGGroup::transform(const MGTransf& tr)//general transformation.
{
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++)
		(**i).transform(tr);

}

//Construct a null newed MGGroup from the type id TID.
MGGroup* MGNullGroup(long TID){
	MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
	return static_cast<MGGroup*>(reg->create_gel(TID));
}

AUTO_GEL_REGISTER(MGGroup, MGGROUP_TID);
AUTO_GEL_REGISTER(MGAppearance, MGAPPEARANCE_TID);

