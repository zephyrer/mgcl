/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPosition_HH_
#define _MGPosition_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/Vector.h"

// MGPosition.h
// Header for MGPosition.

// Forward Declaration.
class MGMatrix;
class MGTransf;
class MGPoint;
class MGCurve;
class MGSurface;
class MGCParam_list;
class MGPosition_list;
class MGIfstream;
class MGOfstream;

/// Represent a positional data.
class MGCLASS MGPosition{

public:

///Translation of the position.
MGDECL inline friend 
MGPosition operator+(const MGPosition& p1,const MGVector& vec){return p1.m_element+vec;};
MGDECL inline friend 
MGPosition operator+(const MGPosition& p1,const MGPosition& p2){return p1.m_element+p2.m_element;};
MGDECL inline friend 
MGPosition operator+ (const MGVector& v, const MGPosition& p){return v+p.m_element;};

/// 自身のPositionと与えられたPositionの減算してMGVectorを生成
MGDECL inline friend 
MGVector operator-(const MGPosition& p1,const MGPosition& p2){return p1.m_element-p2.m_element;}

/// 自身のPositionと与えられたVectorの減算してMGPositionを生成
MGDECL inline friend 
MGPosition operator-(const MGPosition& p1,const MGVector& v){return p1.m_element-v;}

/// 自身のPositionと与えられたVectorの減算してMGPositionを生成
MGDECL inline friend 
MGPosition operator-(const MGVector& v,const MGPosition& p1){return v-p1.m_element;}

/// 自身の点と与えられたベクトルの内積を行う
///Inner product of a osition and a vector.
MGDECL inline friend 
double operator%(const MGPosition& p1,const MGVector& v){return p1.m_element%v;};
	
/// MatrixによるPositionの変換を行いObjectを生成
///Matrix transformation of the position.
MGDECL friend 
MGPosition operator*(const MGPosition& p1,const MGMatrix& mat);

/// PositionのTransformを行いVectorを生成
///General transformation of the position.
MGDECL friend 
MGPosition operator*(const MGPosition& p1,const MGTransf& tr);

/// PositionのScalar乗算を行いObjectを生成
///Scaling of the position.
MGDECL inline friend 
MGPosition operator*(double s, const MGPosition& p){return p*s;};

/// Scalarの乗算を行いPositionを生成
MGDECL inline friend 
MGPosition operator*(const MGPosition& p1,double s){return p1.m_element*s;};

/// Scalar除算を行いObjectを生成
///Scaling of the position.
MGDECL inline friend 
MGPosition operator/(const MGPosition& p1,double s){return p1.m_element/s;};

///Debug Function
MGDECL friend std::ostream& operator<<(std::ostream&, const MGPosition&);

/// 与えられたPositionの成分の値を比較し、同じであれば TRUE を返却
///Comparison of two positions.
MGDECL friend bool operator==(const MGPosition& p1,const MGPosition& p2);
MGDECL friend bool operator==(const MGVector& p1,const MGPosition& p2);
MGDECL friend bool operator==(const MGPosition& p1,const MGVector& p2);
MGDECL inline friend 
bool operator!=(const MGPosition& p1,const MGPosition& p2){return !(p1==p2);};
MGDECL inline friend 
bool operator!=(const MGVector& p1,const MGPosition& p2){return !(p1==p2);};
MGDECL inline friend 
bool operator!=(const MGPosition& p1,const MGVector& p2){return !(p1==p2);};

///Test if this position is less than p2.
///Comparison depends on two positions' length.
MGDECL inline friend 
bool operator<(const MGPosition& p1,const MGPosition& p2){return p1.m_element<p2.m_element;};
MGDECL inline friend 
bool operator<=(const MGPosition& p1,const MGPosition& p2){return p1.m_element<=p2.m_element;};
MGDECL inline friend 
bool operator>(const MGPosition& p1,const MGPosition& p2){return p1.m_element>p2.m_element;};
MGDECL inline friend 
bool operator>=(const MGPosition& p1,const MGPosition& p2){return p1.m_element>=p2.m_element;};

////////// Constructor コンストラクタ ////////////

/// Conversion Constructor from a point
MGPosition(const MGPoint& point);

/// Conversion Constructor from a vector
MGPosition(const MGVector& vec): m_element(vec){;};

///Void constructor void コンストラクタ
explicit MGPosition(size_t sdim=0):m_element(sdim,0.0){;};

///Construct 2D position by providing x,y coordinate data.
MGPosition(double x, double y):m_element(x,y){;};

///Construct 3D position by providing x,y,z coordinate data.
MGPosition(double x, double y, double z):m_element(x,y,z){;};

///Construct 4D position by providing x,y,z,w coordinate data.
MGPosition(double x, double y, double z, double w):m_element(x,y,z,w){;};

/// double の配列でcoordinate valueを指定しPositionを生成する。
///Construct sdim space dimension positional data providing each
///coordinate data through doble array.
///***** This is the fundamental constructor.*****
MGPosition(size_t sdim, const double* v ):m_element(sdim,v){;};

///Construct position by copying old position, changing space dimension and
///ordering of old coordinates.
MGPosition(size_t sdim, const MGPosition& p
			, size_t start1=0, size_t start2=0)
			:m_element(sdim, p, start1, start2){;};

///Construct from std::vector<double>
MGPosition(const std::vector<double>& darrays): m_element(darrays){;};

///Copy constructor
///	MGPosition ( const MGPosition& ); We use default copy constructor.

//////////// Destructor ////////////
//	~MGPosition(); We use default destructor.

//////////// Operator Oveload ////////////

///Assignment
///Update position data by array of double.
MGPosition& operator=(const double* a);

///Return i-th element of the position.
double operator[] (size_t i) const{return m_element.ref(i);}
double operator() (size_t i) const{return m_element.ref(i);}

///Access to i-th element.
double& operator()(size_t i){return m_element(i);};

/// 自身のPositionに与えられたPositionを加算して自身のPositionとする 
///Translation of the position.
MGPosition& operator+= (const MGVector& vec);
MGPosition& operator+= (const MGPosition& pos);

/// 単項マイナス。自身のPositionを反転し、Objectを生成
///Unary minus. Negate all of the elements of the position.
MGPosition operator- () const;

/// 自身のPositionから与えられたVectorを減算し自身のPositionとする
///Translation of the position.
MGPosition& operator-= (const MGVector& vec);

/// Scalarの乗算を行い自身のPositionとする
///Scaling of the position.
MGPosition& operator*=(double scale);

/// MatrixによるPositionの変換を行い自身のPositionとする
///Matrix transformation of the position.
MGPosition& operator*= (const MGMatrix&);

/// PositionのTransformを行いPositionを生成して，
/// 自身のPositionとする
///General transformation of the position.
MGPosition& operator*= (const MGTransf&);

/// Scalar除算を行い自身のPositionとする
///Scaling of the position.
MGPosition& operator/= (double);

//////////// Member Function ////////////

///Let this be the center of the rotation, then compute the angle rotated 
///around the normal from start to end.
///angle(start,end,normal)+angle(end,start,normal)=2*pai always holds.
double angle(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& normal
)const;

///Clear all the element by the value init.
MGPosition& clear(double init=0.0);

///Compute the closest point parameter value of the curve from this point.
///Function's return value is the parameter value of the curve.
double closest(const MGCurve& curve) const;

///Compute the closest point parameter value (u,v)of the surface
///from this point.
MGPosition closest(const MGSurface& surf) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGPosition* clone() const{return new MGPosition(*this);};

///Return the 1st address of the array of the point double data.
const double* data()const{return m_element.data();};
double* data(){return m_element.data();};

///Display list name.
size_t dlist_name()const{return size_t(this);};

///Return the distance of this and P2.
double distance(const MGPosition& P2)const;

/// Generate a Position by interpolating two Position. Input scalar is
/// a ratio t2. When t2 is zero, output position is a copy of the own position.
///Output=(*this)*(1-t2)+vec2*t2.
MGPosition interpolate(double t2, const MGPosition& vec2) const;

///Test if this, P2, and P3 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGPosition& P2,
	const MGPosition& P3
)const{return m_element.is_collinear(P2,P3);};

///Test if this is null.
bool is_null()const{return m_element.is_null();};

/// Positionと原点との距離を求める。
///Return the lenght between the origin(0,0,0) and the position.
double len() const{ return m_element.len();};

/// 点が曲線上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
///Test if the position is on a curve. If on, return the parameter value.
///Even if not on, return the nearest point of the curve.
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGCurve& curve,	///< Curve
	double& t	///< Parameter value of the nearest point on the curve.
) const;

/// 点が曲面上にあるかを調べる。曲面上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
///Test if the position is on a surface. If on, return the parameter value.
///Even if not on, return the nearest point of the surface.
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGSurface& surf,	///< Surface pointer
	MGPosition& uv	///<Parameter value of the nearest point on surface.
) const;

///PD116=Point.
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Return curve's parameter value of this point.
/// If this point is not on the curve, return the nearest point's parameter
/// value on the curve.
double param(const MGCurve& crv) const;

/// Return surface's parameter value of this point.
/// If this point is not on the surface, return the nearest point's parameter
/// value on the surface.
MGPosition param(const MGSurface& srf) const;

///Compute all foot points of the perpendicular line from this point to
///a curve.
/// ポイントから与曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
MGCParam_list perps(
	const MGCurve& crv		///<Curve
)const;

///Compute all foot points of the perpendicular line from this point to
///a surface.
/// ポイントから与曲面へ下ろした垂線の足の，曲面のパラメータ値を
/// すべて求める。
MGPosition_list perps(
	const MGSurface& srf	///<Surface
) const;

/// 自身の点を原点からこの点までのベクトルとして
///ベクトル(v2)に射影したベクトルを求める。
/// v2 が 零ベクトルのとき(*this)が返る。
MGVector project(const MGVector& v2) const{return m_element.project(v2);};

double ref(size_t i) const{return m_element.ref(i);}

///Resize the position, that is , change the space dimension.
///When this is enlarged, the extra space will contain garbages.
void resize(size_t new_sdim){m_element.resize(new_sdim);};

///Get the space dimension
size_t sdim() const {return m_element.sdim();}

///Set this as a null position.
void set_null(){m_element.set_null();}

///Store vec2 data into *this.
///Store length is minimum of len() and vec2.len().
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	size_t i,				///<Displacement of *this.
	const MGVector& vec2,	///<Vector 2.
	size_t j=0			///<Displacement of vec2.
);

///Store vec2 data into *this.
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	size_t i,				///<Displacement of *this.
	const MGVector& vec2,	///<Vector 2.
	size_t j,				///<Displacement of vec2.
	size_t len			///<Length to store 
);

///swap two coordinates.
///swap coordinates (i) and (j).
void swap(size_t i, size_t j){m_element.swap(i,j);};

const MGVector& vector() const{ return m_element;}

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:
	MGVector m_element;		///<Position is a vector.

};

///Test if P1, P2, and P3 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
);

///Compute the angel around the Normal N in radian range[0., 2*pia).
///Here V1=P1-origin, V2=P2-origin.
///Although N is assumed to be parallel to N2=V1*v2, N may not perpendicular
///to v1 and v2, in which case, the projected normal to N2 is used to measure the angle.
///angle(origin,P1,P2,N)+angle(origin,P2,P1,N)=2*pai always holds.
inline double angle(
	const MGPosition& origin,
	const MGPosition& P1,
	const MGPosition& P2,
	const MGVector& N
){
	return origin.angle(P1,P2,N);
};

/** @} */ // end of BASE group
#endif
