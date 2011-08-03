/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBox_HH_
#define _MGBox_HH_
/** @addtogroup BASE
 *  @{
 */

#include <vector>
#include <assert.h>
#include "mg/Interval.h"

// MGBox.h
// Header for class MGBox.
//

//Forward Declaration
class MGVector;
class MGPosition;
class MGMatrix;
class MGTransf;
class MGStraight;
class MGPlane;
class MGIfstream;
class MGOfstream;

///Defines Box of any space dimendion.
/// MGBox expresses n-dimensional space range using MGInterval, which is
/// one dimension range of double values. All of the MGObject's have box.
class MGCLASS MGBox{

public:

///box translation.
MGDECL friend MGBox operator+(const MGVector& v, const MGBox& b);

///Boxを拡大してできるオブジェクトを生成する.
///Generates a box by scaling.
MGDECL friend MGBox operator* (double scale, const MGBox&);

///Print out Debug Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGBox&);

////////Constructor////////

///Copy constructor.
MGBox(const MGBox& box2);

///Void constructor.
///初期化なしでBoxを生成する。
explicit MGBox(size_t dim=0):m_sdim(0), m_range(0){if(dim) get_area(dim);};

///Construct 2D Box, given x and y intervals.
///ｘ，ｙの各インターバルを指定して2D Boxを生成する。
MGBox(const MGInterval& xspan, const MGInterval& yspan);

///Construct 3D Box, given x, y, and z intervals.
///ｘ，ｙ，ｚの各インターバルを指定して3D Boxを生成する。
MGBox(const MGInterval& xspan, const MGInterval& yspan, const MGInterval& zspan);

///Construct Box, given center point and sizes of each coordinates.
///中心点と各辺の幅を指定しBoxを生成する。
///The size's dimension is center.sdim().
MGBox(const MGPosition& center, double* size);

///Construct Box, given center point and a size of each coordinates.
///中心点と幅を指定しBoxを生成する（すべての辺にたいして同一の値）。
///The size is applied to all the coordinates.
MGBox(const MGPosition& center, double size=0.);

///Construct Box, given two points.
///２点からBoxを生成する。
MGBox(const MGPosition&, const MGPosition&);

///Construct Box which contains both input box and a point.
MGBox(const MGBox&, const MGPosition&);

///Construct Box by providing each Interval.
///****This is the fundamental constructor.****
MGBox(size_t dim, const MGInterval*);

///Construct Box by copying old Box, changing space dimension and
///ordering of old coordinates.
///(*this)(start1)=box(start2), and so on.
MGBox(
  size_t dim,		///<space dimension.
  const MGBox& box,	///<original box.
  size_t start1=0,	///<space dimension subscript of the constructed box for start2.
  size_t start2=0	///<space dimension subscript of the input box.
);

////////Destructor///////////
~MGBox();

////////Operator Overload(演算子定義)//////

///Assignment.
MGBox& operator= (const MGBox& box2);

///Return i-th Inteval.
const MGInterval& operator[](size_t i) const {return ref(i);}
const MGInterval& operator()(size_t i) const {return ref(i);}

///Access to i-th Inteval.
MGInterval& operator[](size_t i){assert(i<sdim());return m_range[i];}
MGInterval& operator()(size_t i){assert(i<sdim());return m_range[i];}

///Generate a box that translates the original one.
///Boxを順方向に平行移動してできるオブジェクトを生成する。 
MGBox operator+ (const MGVector& v) const;

///Update the box by translation.
///Boxを順方向に平行移動し自身のBoxとする。
MGBox& operator+= (const MGVector&);

///Generate a box that translates the original one.
///Boxを逆方向に平行移動してできるオブジェクトを生成する。
MGBox operator- (const MGVector&) const;

///Update the box by translation.
///Boxを逆方向に平行移動し自身のBoxとする。
MGBox& operator-= (const MGVector&);

///Generate a box by scaling the box.
///Boxを拡大してできるオブジェクトを生成する。
MGBox operator* (double) const;

///Generate a box by multiplying matrix to the original one.
///与えられたマトリックスでBoxの変換を行いオブジェクトを生成する。
MGBox operator* (const MGMatrix& ) const;

///Generate a box by multiplying transformation to the original one.
///与えられた変換でBoxの８点のトランスフォームを行い、
///それらを囲むBoxのオブジェクトを生成する。
MGBox operator* (const MGTransf& ) const;

///Update the box by scaling.
///自身のBoxを拡大し自身のBoxとする。
MGBox& operator*= ( double );

///update the box by multiplying matrix.
///与えられたマトリックスでBoxの変換を行い自身のBoxとする。
MGBox& operator*= (const MGMatrix& );

///update the box by multiplying transformation.
///与えられた変換でBoxの８点のトランスフォームを行い、
///それらを囲むBoxを自身のBoxとする。
MGBox& operator*= (const MGTransf& );

///Generate a box by scaling the box.
///Boxを縮小してできるオブジェクトを生成する。 
MGBox operator/ ( double ) const;

///Update the box by scaling.
///Boxを縮小し自身のBoxとする。
MGBox& operator/= ( double );

///Generate a box by or operation(minimum box that includes both boxes).
///自身のBoxと与えられたBoxを内包する最小のBoxを生成する。
MGBox operator| ( const MGBox& ) const;

///Update the box by or operation(minimum box that includes both boxes).
///自身のBoxと与えられたBoxを内包する最小のBoxを
///自身のBoxとする。 
MGBox& operator|= ( const MGBox& );

///Generate a box by and operation (minimum box that is common to both boxes).
///自身のBoxと与えられたBoxの共通部分のBoxを生成する。 
MGBox operator& ( const MGBox& ) const;

///Uppdate the box by and operation (minimum box that is common to both boxes).
///自身のBoxと与えられたBoxの共通部分のBoxを
///自身のBoxとする。
MGBox& operator &= ( const MGBox& );

///Equal and not equal operations of two boxes.
///自身と与えられたBoxが等しいかどうかを返却する。
///一方のBoxの全ての辺がもう一方の相当する辺に重なる場合 
///TRUE を返却する。
bool operator== (const MGBox& box2) const;
bool operator!= (const MGBox& box2) const{return !(*this==box2);}

///与えられたポジションが自身のBox内に含まれているか返却する。
///ポジションがBox内にある場合 TRUE を返却する。
///完全にポジションが自身のBox内にある場合、含まれるとみなす。
///Returns true if the box includes the position.
bool operator >> (const MGPosition&) const;

///Returns true if the box includes the second box.
///自身のBoxが与えられたBoxを囲んでいるか返却する。
///与えられたBoxが自身のBox内にある場合 TRUE を返却する。
///与えられたBoxがNULL の場合は、FALSE 。
bool operator>> (const MGBox&) const;

///Returns true if the box does not includes the position.
///与えられたポジションが自身のBoxの範囲外にあるかどうか返却する。
bool operator<< (const MGPosition& pt) const{return !((*this)>>pt);}

///Returns true if the box does not include the second box.
///与えられたBoxが自身のBoxを囲んでいるかどうか返却する。 
bool operator<< (const MGBox& box2) const{return box2>>(*this);}

////////Member Function////////

///Test if the straight line sl is crossing this box or not.
///Function's return value is true if a part of sl is included in this box,
///false if not.
bool crossing(const MGStraight& sl)const;

///Test if the plane is cutting this box or not.
///Function's return value is true: if the plane is cutting the box,
///false if all of the vertices of the box are in one side of the plane.
bool cutting(const MGPlane& plane)const;

///Compute the distance from a point to the box.
double distance(const MGPosition& P) const;

///Dump Function.
int dump(MGOfstream& ) const;

///Calculate dump size.
size_t dump_size() const;

///Expand the box by MGTolerance::rc_zero().
///That is,
///operator*=(1.+ MGTolerance::rc_zero() ) will be executed.
void expand();

///Expand the box by len.
///This box will be updated so that the center of
///the box will not be moved and the box is widened by len for each coordinate.
///That is, ref(i).high() is set to ref(i).high()+len
///and ref(i).low() is set to ref(i).low()-len for each i.
void expand(double len);

///Expand the box by len[].
///This box will be updated so that the center of
///the box will not be moved and the box is widened by len[i]
///for each coordinate. That is, ref(i).high() is set to ref(i).high()+len[i]
///and ref(i).low() is set to ref(i).low()-len[i] for each i.
void expand(double* len);

///Expand the box so that this contains the position P.
void expand(const MGPosition& P);

///Return ture if Box is empty.
///Boxが empty かどうか返却する
bool empty() const;

///Return true if box is finite.
bool finite() const;

///Return maximum(high), minimum(low), or middle(mid) points of
///the Box.
///Boxの対角線の両端と中心を返却する。全ての座標値が
///最小の点が low () で、最大の点が high () で返却される。 
MGPosition high() const;
MGPosition low() const;
MGPosition mid() const;

///Test if this includes the origin(0., 0., ...).
///原点が自身のBox内に含まれているか返却する。
bool includes_origin()const;

///Test if the point P is included in this box.
bool includes(const MGPosition& P)const{return operator>>(P);};

///Test if this is null box.
bool is_null()const{return m_sdim==0;}

///Return diagonal line length.
///ボックスの対角線長さを求める
double len() const;
double length()const{return len();};

///Reference to i-th Interval.
const MGInterval& ref(size_t i) const;

///Restore Function.
int restore(MGIfstream& );

///Return space dimension.
size_t sdim() const{return m_sdim;};

///Update maximum coordinate values by input data.
///自身のBoxの最大座標値を指定された点に変更する.
MGBox& set_high(const MGPosition&);

///Update minimum coordinate values by input data.
///自身のBoxの最小座標値を指定された点に変更する.
MGBox& set_low(const MGPosition&);

///Set this box as a null box.
void set_null();
														  
///Return the type of i-th interval.
MGINTERVAL_TYPE type(size_t i) const{
	if(i>=m_sdim) return MGINTERVAL_EMPTY;
	return m_range[i].type();
}

/// vertex computes all the vertices of the box.
/// Return array of Position as the vertices.
std::vector<MGPosition> vertex() const;

private:
///Member Data
	MGInterval m_rData[3];	///<Coordinate Range work area when m_sdim<=3.
	size_t m_sdim;			///<Space dimension of this box,
			///<m_range's length may differ from m_sdim,
			///<That is, m_sdim<=m_range's length.
	MGInterval* m_range;///<If m_sdim<=3, the area of m_rData will be used,
						///<If m_sdim>3, newed area will be used.

///Test if the straight line sl is crossing this box or not.
///Function's return value is true if a part of sl is included in this box,
///false if not. 2D version of crossing().
bool crossing2D(const MGStraight& sl)const;

///Get the area of m_range for the space dimension sdim.
///Result area will contain garbages.
///get_area will use m_sdim as input. m_sdim must be valid.
void get_area(size_t dim);

///Resize the m_range.
///If sim> current sdim(), the old data are guaranteed to hold in resized box.
void resize(size_t dim);

};

/** @} */ // end of BASE group
#endif
