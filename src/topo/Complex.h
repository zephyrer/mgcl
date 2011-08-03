/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGComplex_HH_
#define _MGComplex_HH_

#include "mg/Box.h"
#include "topo/Topology.h"

#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<MGCellNB*>;
#pragma warning( pop )
#else
#include <list>
#endif

//
//Define MGComplex Class.

class MGPosition;
class MGCellNB;
class MGCellMap;
class MGBoundary;

/** @addtogroup TOPO
 *  @{
 */

///MGComplex is a container of parameter cells and binder cells.
class MGCLASS MGComplex: public MGTopology{

public:
#if defined(MGCL_DLL)
typedef MGListProxy<MGCellNB*> container_type;
#else
typedef std::list<MGCellNB*> container_type;
#endif

typedef container_type::iterator cellItr;
typedef container_type::const_iterator const_cellItr;

typedef container_type::iterator pcellItr;
typedef container_type::const_iterator const_pcellItr;

typedef container_type::iterator bcellItr;
typedef container_type::const_iterator const_bcellItr;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGDECL friend MGComplex operator* (double s, const MGComplex& complex);

/////////Constructor/////////

///Void constructor.
MGComplex();

///Copy constructor. Copy as a boundary complex of parent.
///When parent is not specified, this is ordinary world complex.
/// Not a boundary complex.
MGComplex(const MGComplex& complex);

///Construct of one cell.
///The second form takes the ownership of the cell, must be newed object.
MGComplex(const MGCellNB& cell);
MGComplex(MGCellNB* cell);

/////////Destructor/////////

virtual ~MGComplex();

/////////operator overload/////////

///Assignment.
///When the leaf object of this and topo2 are not equal, this assignment
///does nothing.
virtual MGComplex& operator=(const MGGel& gel2);
virtual MGComplex& operator=(const MGComplex& gel2);

/// Complexに平行移動を行ないオブジェクトを生成する。
///Translation of the Complex
MGComplex operator+ (const MGVector& v) const;

/// Complexに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Complex
MGComplex operator- (const MGVector& v) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGComplex operator* (double s) const;

/// 与えられた変換でComplexの変換を行い，Complexを作成する。
///Transformation of the Complex by a matrix.
MGComplex operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないComplexを生成する。
///Transformation of the Complex by a MGTransf.
MGComplex operator* (const MGTransf& tr) const;

///Object transformation.
virtual MGComplex& operator+=(const MGVector& v);
virtual MGComplex& operator-=(const MGVector& v);
virtual MGComplex& operator*=(double scale);
virtual MGComplex& operator*=(const MGMatrix& mat);
virtual MGComplex& operator*=(const MGTransf& tr);

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGComplex operator/ (double s) const{return (*this)*(1./s);};

///comparison
virtual bool operator<(const MGComplex& gel2)const;
virtual bool operator<(const MGGel& gel2)const;

/////////Member Function/////////

///Obtain first bcell iterator.
const_bcellItr bcell_begin() const{ return m_bcells.begin();};
bcellItr bcell_begin() { return m_bcells.begin();};

///Obtain end bcell iterator(next of the last bcell).
const_bcellItr bcell_end() const{ return m_bcells.end();};
bcellItr bcell_end() { return m_bcells.end();};

///Obtain i-the pcell in the m_bcells sequence. MGCellNB version.
const MGCellNB* bcelli(size_t i) const;
MGCellNB* bcelli(size_t i);

///Obtain i-the pcell in the m_bcells sequence. MGCellNB version.
const_bcellItr bcellIterator(size_t i) const;
bcellItr bcellIterator(size_t i);

///Cehck if bcell exist.
bool bcell_exist()const;

///Obtain the binder of i-the pcell, may be null.
MGCellNB* binder(size_t i) const;

///Obtain binders of all the pcells of the boundary.
///i-th binder of the fucntion's returned value is the binder
///of i-th pcell of the boundary, and may be null.
std::vector<MGCellNB*> binders() const;

///Return the box of this complex.
const MGBox& box() const;

///Compute barycenter of all the vertex(binder cell of 0D manifold
///dimension).
MGPosition center() const;

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGComplex* clone() const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D curve in the topology's star cell world coordinates.
///The object is converted to curve(s) and is drawn.
virtual void drawWire_in_star(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void draw3DVertex()const;

///Draw 3D point(vertex) in star cell's world coordinates.
///The object is converted to point(s) and is drawn.
virtual void draw3DVertex_in_star()const;

///Erase first pcell element, including binder cells
///that will be free when the first pcell is erased.
void erase_first_pcell();	///erase first pcell.

///Erase last pcell element, including binder cells
///that will be free when the last pcell is erased.
void erase_last_pcell();	///erase last pcell.

///Erase the pcell element, including binder cells
///that will be free when the last pcell is erased.
void erase_pcell(MGCellNB* pcell);

///Get fisrt pcell pointer.
const MGCellNB* first_pcell()const;
MGCellNB* first_pcell();

///Return Object's type ID (TID)
virtual long identify_type()const;

///Test if this complex includes the MGCellNB cell as a contituent.
///Returns true if cell is included in this complex.
bool includes(const MGCellNB* cell)const;

///Get last pcell pointer.
const MGCellNB* last_pcell()const;
MGCellNB* last_pcell();

///Get manifold dimension.
virtual unsigned manifold_dimension() const;

///count number of bcells of the complex
size_t number_of_bcells() const{return m_bcells.size();};

///count number of pcells of the complex
size_t number_of_pcells() const{return m_pcells.size();};

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Obtain first pcell iterator.
const_pcellItr pcell_begin() const{ return m_pcells.begin();};
pcellItr pcell_begin() { return m_pcells.begin();};

///Obtain end pcell iterator(next of the last pcell).
const_pcellItr pcell_end() const{ return m_pcells.end();};
pcellItr pcell_end() { return m_pcells.end();};

///Obtain i-the pcell. MGCellNB version.
const MGCellNB* pcelli(size_t i) const;
MGCellNB* pcelli(size_t i);

///Obtain i-th pcell iterator.
const_pcellItr pcellIterator(size_t i) const;
pcellItr pcellIterator(size_t i);

///Cehck if pcell exist.
bool pcell_exist()const;

///Obtain pcells that constitute the boundary.
///Let pcellvec[.] be pcells' return value and bindervec[.] be
///binders' return value. Then pcellvec[i] corresponds to bindervec[i].
///bindervec[i] is binder cell of i-th pcell element of the boundary.
std::vector<MGCellNB*> pcells();
std::vector<const MGCellNB*> pcells() const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Get star cell pointer if this complex is a boundary of a cell.
///Else, null will be returned.
const MGCellNB* star() const;
MGCellNB* star();

protected:
	mutable MGBox m_box;	///<Box of the complex.
		///<Initially this is null, and will be computed if necessary.
		///<***Currently Box is computed only from pcells, not correct box.

///Fundamental constructor.
///Construct from a list of pcells.
///This constructor takes the ownership of all pcells in pcells.
///Fundamental constructor of MGComplex does not take box as input since
///this constructor has to do binder append treatment that is attached to
///pcells.
explicit MGComplex(std::list<MGCellNB*>& pcells);

///Binder cells of the pcells in complex will be registered in cmap.
MGComplex(const MGComplex& complex,
	  	MGCellMap& cmap);		///cellmap to register binder association.

///Append a binder to the end of bcell sequence.
cellItr append_bcell(MGCellNB* cell)const;

///Insert a Cell(binder or parameter) before the position loc.
///loc is the iterator of m_pcells or m_bcells according to pcell.
cellItr add_cell(
	MGCellNB* cell,	///<Cell to insert
	bool pcell,		///<indicates if input cell is parameter or binder cell.
	cellItr loc	///<Iterator that indicates insert position. Insert cell befoer 
					///<the iterator loc.
);

///Append a PCell to the end of pcell sequence.
cellItr append_pcell(MGCellNB* cell);

///Compute the box from the scratch.
virtual void compute_box() const;

///Copy all pcells of comp into this,
///but does not copy binders of the pcells.
void copy_without_binders(const MGComplex& comp);

///Erase all elements(PCells and BCells).
///erase_all_elements does free and destruct all the elements.
///erase_all_elements does not maintain m_box.
void erase_all_elements();

///Prepend a PCell to the end of pcell sequence.
cellItr prepend_pcell(MGCellNB* cell);

virtual std::string whoami()const{return "Complex";};

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///Assighment.
///When the leaf object of this and topo2 are not equal, this assignment
///does nothing.
MGComplex& set_complex(const MGComplex& comp2);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

private:

	container_type m_pcells;	///<list of pcell elements.
	mutable container_type m_bcells;///<vector of bcell elements.
	
///****All of the following transformation operators are applied to
///binders of the boundary's pcells, not to pcell themselves.****

/// BoundaryのBinderに平行移動を行ない自身のBoundaryとする。
///Translation of the boundary.
void binderTr(const MGVector& v);

///BoundaryのBinderにスケーリングを行い自身のBoundaryとする。
///Scaling of the boundary by a double.
void binderTr(double s);

/// 与えられた変換でBoundaryのBinderに変換を行い自身のBoundaryとする。
///Transformation of the boundary by a matrix.
void binderTr(const MGMatrix& mat);

/// 与えられた変換によってBinderにトランスフォームをおこない自身のBoundaryにする。
///Transformation of the boundary by a MGTransf.
void binderTr(const MGTransf& tr);

///Copy all the pcells of comp into this.
void copy_all_elements(const MGComplex& comp);

///Copy all the pcells of comp into this.
///Binder cells of the pcells in comp will be registered in cmap.
void copy_all_elements(
	const MGComplex& comp,
	MGCellMap& cmap		///<cellmap to register binder association.
);

///Erase bcell element if it is registered in bcell sequence of this complex.
///Destruct the bcell if found in this complex.
void erase_binder(const MGCellNB*);

friend class MGCellBase;
friend class MGCellNB;
friend class MGCell;
friend class MGEdge;
};

/** @} */ // end of TOPO group
#endif
