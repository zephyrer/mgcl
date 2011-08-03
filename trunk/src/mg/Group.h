/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGroup_HH_
#define _MGGroup_HH_

#include "mg/Plist.h"
#include "mg/Pvector.h"
#include "mg/Object.h"
#include "mg/Gel.h"
#include "mg/AttribedGel.h"

//
//Define MGGroup Class.
//
class MGIfstream;
class MGOfstream;
class MGBox;
class MGObject;
class MGAttrib;
class MGGLAttrib;
class MGAppearance;
class MGContext;

/** @addtogroup GelRelated
 *  @{
 */

///MGGroup is a class which constains MGGel elements.
///MGGroup provides functions:
///(1) a container of MGGel(MGAtrribedGel(MGObject or MGGroup), MGAttrib) elements.
///(2) To attach an appearance to the context of the group.
class MGCLASS MGGroup:public MGAttribedGel{

public:
///	別名定義
	typedef MGPlist<MGGel>::iterator iterator ;
	typedef MGPlist<MGGel>::const_iterator const_iterator;
	typedef MGPlist<MGGel>::reverse_iterator reverse_iterator ;
	typedef MGPlist<MGGel>::const_reverse_iterator const_reverse_iterator;

///Public member data. A list of MGGel*.
		MGPlist<MGGel> m_gels;

////////////Constructor////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGGroup();
MGGroup(MGGel* gel);

///Construct MGGroup from the file made by MGOfstream.
///error contains error flag as:
///=0: open succeeded.
///=1: file not found, or could not be opened.
///=2: file found, but, the format is not MGCL format.
///When error!=0, MGGroup contains no MGGel's.
MGGroup(const char* file, int& error);

///Copy constructor.
///MGGroup(const MGGroup& obj2);

///Virtual Destructor
///virtual ~MGGroup();

//////////////////operator overloaded////////////////

///Assignment.
///When the leaf object of this and gel2 are not equal, this assignment
///does nothing.
virtual MGGroup& operator=(const MGGel& gel2);
virtual MGGroup& operator=(const MGGroup& gel2);

///comparison
virtual bool operator<(const MGGroup& gel2)const;
virtual bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Get the MGAppearance pointer in this group. If not defined, null will be
///returned.
MGAppearance* appearance();
const MGAppearance* appearance()const;

/// Return(but does not remove) last element in the group.
/// If list is empty, behavior is undefined.
const MGGel* back() const{return m_gels.back();};
MGGel* back(){return m_gels.back();};

/// Return const_iterator at the beginning of list.
const_iterator begin() const{return m_gels.begin();}
iterator begin(){return m_gels.begin();}

///Get the box of the group.
///If no objects were included in this group, null box will be returned.
MGBox box()const;

/// clear list, that is, erase all the elements in the MGGel.
void clear(){m_gels.clear();};

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
virtual MGGroup* clone()const;

///Get the MGContext pointer stored in this group. If not defined, null will be
///returned.
MGContext* context();
const MGContext* context()const;

///Delete a display list of this gel.
void delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	///<mgSysGL pointer will be apppended.
)const;

///////display member function.
virtual void display_arrows()const;
virtual void display_break_points()const;
virtual void display_control_polygon()const;
virtual void display_curvatures(
	double	scale,	///scaling of the graph.
	int		density,///densitiy of the graph.
	bool	use_radius///true:radius display, false:curvature display.
)const;

///make this group has appearance and get the MGAppearance pointer.
MGAppearance* ensure_appearance();

///Return true (1) if there are no items in the MGGel,
/// false(0) otherwise.
bool empty() const{return m_gels.empty();};

/// Return const_iterator at the end of MGGel.
const_iterator end() const{return m_gels.end();};
iterator end(){return m_gels.end();};

/// erase element x. Function's return value is the following iterator
/// of the erased element x.
iterator erase(iterator x){return m_gels.erase(x);};

/// erase sequence [first, last). Function's return value is the following iterator
/// of the erased elements.
iterator erase(iterator first, iterator last){return m_gels.erase(first,last);};

///Find the position of the gel in the gel list and the group pointer
///which includes the gel. Searching will be done into the member group gel
///of this list.
const_iterator find(
	const MGGel* gel,
	const MGGroup*& grp	///<The group that includes gel will be returned.
						///<When not found grp=null will be returned
)const;
iterator find(
	MGGel* gel,
	MGGroup*& grp	///<The group that includes gel will be returned,
					///<When not found grp=null will be returned.
);

/// Return(but does not remove) first element in the MGGel.
/// If list is empty, behavior is undefined.
const MGGel* front() const{return m_gels.front();};
MGGel* front(){return m_gels.front();};

///Return MGGroup pointer if this MGGel is an MGGroup, else return null.
MGGroup* group(){return this;};
const MGGroup* group()const{return this;};

/// Return This object's typeID
virtual long identify_type() const{return MGGROUP_TID;};

///Test if this gel includes an object.
///Function's return value is the 1st object found in the gel list of this
///if this includes an object. Otherwise null will be returned.
const MGObject* includes_object()const;
MGObject* includes_object();

///insert an element x before the position it.
///Function's return value is the iterator of x after inserted.
iterator insert(iterator it, MGGel* x){return m_gels.insert(it,x);};

///Make an MGOfstream file which contains this group.
///The file generated by make_file can be retrieved by the constructor
///MGGroup(const char* file, int error);
///Function's return value is:
///=0: the file is successfully made.
///=1: file could not be opened.
int make_file(const char* file);

///Make a display list of this gel.
///Return is the display list name.
size_t make_display_list(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Make a display list without color of this gel.
///Return is the display list name.
size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Make a display list of only call lists of this group.
///Return is the display list name.
///When gels in this group are included in gels_to_delete, the display list
///will not be made.
size_t make_only_call_list(
	bool no_color=false,///<if true, color attribute will be neglected.
	const std::vector<const MGGel*>* gels_to_delete=0
)const;

///Returns the size of maximum size.
size_t max_size() const{return m_gels.max_size();};

///Get the number of objects included in thie group.
int num_of_objects() const;

///IGES output function.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

/// pop last element.
void pop_back(){m_gels.pop_back();};

/// pop first element.
void pop_front(){m_gels.pop_front();};

void push_appearance(MGAppearance* appr);
void push_context(MGContext* cntx);

/// push element x at the end.
void push_back(MGGel* x);

///push elements in objs at the end. All of the object pointers are
///transfered to this. On return, objs will have no object pointer in it.
void push_back(MGPvector<MGObject>& objs);
void push_back(MGPlist<MGObject>& objs);
void push_back(MGPvector<MGGel>& gels);
void push_back(MGPlist<MGGel>& gels);

/// push element x at the first.
///x must be a newed object, and the ownership will be transfered to thisMGGroup.
void push_front(MGGel* x);

/// Return const_reverse_iterator at the beginning of list.
const_reverse_iterator rbegin() const{return m_gels.rbegin();};
reverse_iterator rbegin(){return m_gels.rbegin();};

/// Return const_reverse_iterator at the end of list.
const_reverse_iterator rend() const{return m_gels.rend();};
reverse_iterator rend(){return m_gels.rend();};

///Release the gel at the position i.
///Returned will be the position after the relesed gel.
iterator release(iterator i){return m_gels.release(i);};

///Remove the MGAppearance of this MGAttribedGel.
void remove_appearance();

///Remove the MGGel* and return the MGGel*. If i is not valid, 
/// behavior is undefined.
///現在の所有権を放棄し、ただのポインタを返す。
MGGel* removeAt(iterator x){return m_gels.removeAt(x);};
	
///reverse sequence.
void reverse(){m_gels.reverse();};

/// Return the number of items that are in the list.
size_t size() const{return m_gels.size();};

///Transform the gel by the argument.

///translation
virtual void transform(const MGVector& v);

///scaling.
virtual void transform(double scale);

///matrix transformation.
virtual void transform(const MGMatrix& mat);

///general transformation.
virtual void transform(const MGTransf& tr);

virtual std::string whoami()const{return "Group";};

protected:

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

private:

///Make a display list of only call lists of this group by the name to push for pick
///and display list name glname.
///When gels in this group are included in gels_to_delete, the display list
///will not be made.
void make_only_call_list_sub(
	size_t name,	///<the name to glPushName for pick
	int mode,	///<=0 when wire, =1 when wire and no color, =2 when SHADING.
	const std::vector<const MGGel*>* gels_to_delete
)const;

friend class MGIfstream;
friend class MGOfstream;
};

///Construct a null newed MGGroup from the type id TID.
MGDECL MGGroup* MGNullGroup(long TID);

/** @} */ // end of GelRelated group
#endif
