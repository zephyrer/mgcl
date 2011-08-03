/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGRLBRep_HH_
#define _MGRLBRep_HH_

#include "mg/LBRep.h"

// MGRLBRep.h
//

// Forward Declaration
class  MGPosition;
class  MGKnotArray;
class  MGCParam_list;
class  MGPosition_list;
class  MGIfstream;
class  MGOfstream;

/** @addtogroup GEO
 *  @{
 */

/// Defines Rational Line B-Representation.
/// This NURBS is a homogeneous form, i.e., B-Coefficients have
/// weight included values. 
/// When usual NURBS form is (xi, yi, zi, wi) ,
/// MGRLBRep form is (xi*wi, yi*wi, zi*wi, wi) for i=0,..., n-1.
class MGCLASS MGRLBRep: public MGCurve {

public:
///Friend Function
MGDECL friend MGRLBRep operator+ (const MGVector& v, const MGRLBRep& lb);
MGDECL friend MGRLBRep operator* (double scale, const MGRLBRep&);

//////////// Constructor ////////////

///Default(dummy) constructor.
MGRLBRep():MGCurve(){;}

///Construct Line NURBS, providing all the member data.
///***** This is the fundamental constructor(when homogeneous=1).*****
MGRLBRep(
	const MGKnotVector& t,		///<Knot Vector.
	const MGBPointSeq& bcoef,	///<Line B-Coef, each of coefficients
		///<includes weight multiplied when homogeneous=true(1),
		///<and not includes when homogeneous =false.
		///<Mximum space dimension id of bcoef is for weight of the rational.
	int homogeneous=1);

///Construct Line NURBS, providing all the member data.
MGRLBRep(
	const MGKnotVector& t,	///<Knot Vector.
	const MGBPointSeq& bcoef,///<Line B-Coef, each of coefficients does not include weights data.
	const std::vector<double>& weights
);

/// Construct ellipse NURBS.
explicit MGRLBRep(const MGEllipse& ellipse);///Original ellipse.

/// Construct a conic section NURBS.
///This conic is defined by ths start and end point, and each tangent,
///and mid-point of the conic.
MGRLBRep(
	const MGPosition& P0, const MGVector& T0,
							///<Start point and its tangent
	const MGPosition& P,	///<Mid point of the conic section
	const MGPosition& P2, const MGVector& T2
							///<End point and its tangent
);

///Approximate an original NURBS by a new knot configuration.
///The new knot config must be inside the range of the original NURBS
///parameter. However new knots may be coarse or fine.
MGRLBRep(
	const MGRLBRep& old_brep,///<Original NURBS.
	const MGKnotVector& t,	///<knot vector
	int &error			///<Error flag.
);

///**** Conversion Constructor.****

/// Convert from Non ratoinal form to Rational form.
///When homogeneous==true(non zero), brep is homogeneous form MGLBRep.
///When homogeneous==false(zero), brep is ordinary MGLBRep and
///will be converted to MGRLBRep. That is, weight=1 elements will be
///added as last space dimension element.
///***** This is the fundamental constructor. *****
MGRLBRep(
const MGLBRep& brep,	///<Original LBRep. This can be ordinary LBRep, or 
	///<homogeneous form of MGRLBRep. When homogeneous form,
	///<the last space dimension elements are weights.
int homogeneous=0		///<true(non zero): homogeneous form,
						///<false(zero):ordinary LBRep.
);

///Gets new NURBS by adding knots to an original NURBS.
MGRLBRep(
	const MGRLBRep& old_brep,	///<Original NURBS.
	const MGKnotArray& knots	///<Knots to add.
);

///Gets new NURBS by connecting two NURBS to one.
MGRLBRep(
	const MGRLBRep& brep1,	///<NURBS 1.
	int continuity,			///<continuity.
	int which,	///<which point of brep1 to which of brep2,
				///<meaingfull when continuity>=0,
			///< =0: start of this and start of brep1,
			///< =1: start of this and end of brep1,
			///< =2: end of this and start of brep1,
			///< =3: end of this and end of brep1,
			///< continuity and which can be obtained using continuity().
	const MGRLBRep& brep2	///NURBS 2.
);

/// Gets new NURBS by computing a part of the original. New one is exactly
/// the same as the original except that it is partial.
///If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim() and k=order().
///Both t1 and t2 must be inside te range of old_brep.
MGRLBRep(
	double t1, double t2, ///<New parameter range. t1 must be less than t2.
	const MGRLBRep& old_brep,	///<Original NURBS.
	int multiple=0 ///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

/// Construct a Line NURBS by changing space dimension and ordering of
///coordinates.
MGRLBRep(
	size_t dim,				///< New space dimension.
	const MGRLBRep& lbrep,	///< Original Line B-rep.
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of original line.
);

///Construct 2D ellipse RLBRep, whose center is origin.
///The ellipse is expressed as below using parameter t.
/// x(t)=a*cos(t),  y(t)=b*sin(t),   angle1<=t<=angle2
MGRLBRep(
	double a, double b,
	double angle1, double angle2
);

///	MGRLBRep(const MGRLBRep&);  ///Copy constructor.
///  We can use default copy constructor.

//////////// Destructor /////////
///	~MGRLBRep();	///We can use default destructor.

//////////// Operator overload ////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGRLBRep& operator=(const MGGel& gel2);
MGRLBRep& operator=(const MGRLBRep& gel2);

///Transformation object construction
MGRLBRep operator+ (const MGVector& ) const;
MGRLBRep operator- (const MGVector& ) const;
MGRLBRep operator* (double) const;
MGRLBRep operator* (const MGMatrix& ) const;
MGRLBRep operator* (const MGTransf& ) const;

///Object transformation.
MGRLBRep& operator+=(const MGVector& v);
MGRLBRep& operator-=(const MGVector& v);
MGRLBRep& operator*=(double scale);
MGRLBRep& operator*=(const MGMatrix& mat);
MGRLBRep& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGRLBRep& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGRLBRep& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator==(const MGLBRep& gel2)const;

//////////// Member Function ////////////

///Returns NURBS Dimension.
size_t bdim() const{return m_line.bdim();}
	
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
///Return minimum box that includes the partial line.
MGBox box_limitted(const MGInterval& l) const;

///Return minimum box that includes the whole line.
MGBox box_unlimit() const;

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,			///<Parameter value for the start of original. 
	double t2			///<Parameter value for the end of original. 
){
	m_line.change_range(t1,t2);
	update_mark();
};

///Changing this object's space dimension.
MGRLBRep& change_dimension(
	size_t dim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Change order of the NURBS. When new order is greater than the original,
///new B-rep is guaranteed to be the same line as the original. However,
///if new order is less than the original one, new line is not the same
///in general.
MGRLBRep& change_order(
unsigned order		///<New order number. 
){m_line.change_order(order); return *this;};

///Access to (i,j)th element of coef
///( left-hand side version)
double& coef(size_t i, size_t j){update_mark(); return m_line.coef(i,j);};		

///Access to (i,j)th element of coef
///(right hand side version)
double coef(size_t i, size_t j) const{return m_line.coef(i,j);};

///Extract (i,j)elements for 0<=j<=sdim() as a vector
///of (sdim()+1) space dimension. The last elemnt is weight.
MGVector coef(size_t i) const{return m_line.coef(i);};

///Returns a pointer to the line b-coef data.
const double* coef_data(size_t i=0, size_t j=0) const
{return m_line.coef_data(i,j);};

///Connect brep2 to this brep to make one B-Representation.
///This parameter range will not be changed, instead brep2's range
///will be so changed that brep2 has the same 1st derivative magnitude
///as the original this brep's at the connecting point(start or end point of
///this).
///continuity and which can be obtained using the fucntion continuity().
void connect(
	int continuity,	///<continuity. must be>=0.
	int which,		///<which point of this to which of brep2.
				///< =0: start of this and start of brep2.
				///< =1: start of this and end of brep2.
				///< =2: end of this and start of brep2.
				///< =3: end of this and end of brep2.
	const MGRLBRep& brep2	///<B-Rep 2.
);

///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two lines are discontinuous.
///  0: G0 continuity, i.e. two lines are connected,
///     but tangents are discontinuous
///  1: C1 continuity, i.e. 1st deriv's  are also continuous,
///     when weights are so arranged.
///  2: C2 continuity, i.e. 2nd deriv's  are also continuous,
///     when weights are so arranged.
int continuity(
	const MGRLBRep& brep2,
	int& which,	///<Indicates which point of this is connected
				///< to which of brep2, is meaingfull when continuity>=0,
				///< =0: start of this to start of brep2,
				///< =1: start of this to end of brep2,
				///< =2: end of this to start of brep2,
				///< =3: end of this to end of brep2.
	double& ratio ///< Ratio of 1st derivatives of the two line will
			///< be returned,
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGRLBRep& coordinate_exchange(size_t i, size_t j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGRLBRep* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const{return clone();};

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGRLBRep* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is newed object and must be deleted.
MGCurve* copy_limitted(const MGInterval& prange) const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	MGPvector<MGCurve>& crv_list,	///<divided curves will be appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
) const;

///Display control polygons using mgGDL::MGDrawPointSeq()
void display_control_polygon()const;

///Draw this line's 1st and 2nd coordinates in 2D space
///using drawing function moveto( , ) and lineto( , ).
///wind[] is the window of the screen to draw the line in.
///Clipping will be performed about the wind[].
///(wind[0], wind[1]) is the center coordinates of the window.
///wind[2] is width and wind[3] is hight of the window. When wind[2]<=0,
///no clipping is performed. Even when wind[2]<=0, wind[3] is necessary 
///to input to specify the resolution of the line. In this case,
///wind[0] and wind[1] are not referended.
///ynum is the resolution of the line, is the number of
///straight line segments for the curve length of wind[3](height of window).
///***draw_2D does not perform box including judment, always performs clipping
///operation and draws the line. Users must do obvious box inclusion test
///if maximum drawing performance is necessary.
void draw_2D(void (*moveto)(int, int), void (*lineto)(int, int),
	const double wind[4],	///<window box to draw the line in.
	size_t ynum			///<Resolution of the line.
)const;

void draw_2D(void (*moveto)(float, float), void (*lineto)(float, float),
	const double wind[4],	///<window box to draw the line in.
	size_t ynum			///<Resolution of the line.
)const;

void draw_2D(void (*moveto)(double, double), void (*lineto)(double, double),
	const double wind[4],	///<window box to draw the line in.
	size_t ynum			///<Resolution of the line.
)const;

///Draw this line's coordinate'th coordinate in 2D space as
///(t, LBRep(coordinate)) when t_is_x is true, 
///or as ( LBRep(coordinate),t)  when t_is_x is false,  
///Here t is the parameter of the LBRep.
///using drawing function moveto(int, int) and lineto(int,int).
///The other behaviours are the same as draw_2D.
void draw_1D(void (*moveto)(int, int), void (*lineto)(int, int),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_1D(void (*moveto)(float, float), void (*lineto)(float, float),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_1D(void (*moveto)(double, double), void (*lineto)(double, double),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum	///<Resolution of the line.
)const;

void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve.
						///<Draw will be performed from t0 to t1.
)const;

/// Evaluate right continuous n'th derivative data.
/// nderiv=0 means positional data evaluation.
MGVector eval(
	double t,		///< Parameter value.
	size_t nderiv=0,///< Order of Derivative.
	int left=0		///<Left continuous(left=true)
					///<or right continuous(left=false).
) const;

///Compute position, 1st and 2nd derivatives.
/// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
void eval_all(
	double tau,			///<Input parameter value(パラメータ値)
	MGPosition& P,		///<Position(位置)
	MGVector& V1,		///<1st derivative(1次微分値)
	MGVector& V2		///<2nd derivative(2次微分値)
) const;

///Evaluate all of i'th derivative data for 0<=i<=nderiv.
///Output will be put on deriv[j+i*sdim()]
///for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
///deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void eval_all(
	double tau,		///< Parameter value to evaluate.
	size_t nderiv,	///< Order of Derivative.
	double* deriv,	///< Output area of size (nderiv+1)*sdim().
	int left=0		///<Left continuous(left=true)
					///<or right continuous(left=false).
) const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Modify the original NURBS by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
MGRLBRep& extend(
	int start,			///<Flag of start or end poit of the line,
						///<If start is true extend on the start point.
	double length,		///<chord length to extend. 
	double dk=0.         ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.,
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point,
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Extrapolate the curve by the parameter value.
MGRLBRep& extend_with_parameter(
	double tau,	///<The parameter value at the end of extended point,
				///<When tau<param_s(), extension will be done at the starting point,
				///<When tau>param_e(), extension will be done at the end point.
	double dk   ///<Coefficient of how curvature should vary at the connecting point.
				///<See extend();
);

///Extracts control points.
///Fucntion's return value is 
///true if control points was obtained, false if not.
bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const;

///Return homogeneous Line B-Representation of the rational B-Spline.
const MGLBRep& homogeneous() const {return m_line;}
MGLBRep& homogeneous(){update_mark();return m_line;}

/// Return This object's typeID
long identify_type() const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum() const;

///Test if this cure is co-planar with the 2nd curve curve2.
///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

/// Spline と Curve の交点を求める。
///Intersection point of spline and curve.
MGCCisect_list isect(const MGCurve&) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;

///Intersection of Spline and Surface

///Intersection of MGRLBRep and Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Compute intersection point of 2D sub NURBS of original B-rep.
///Parameter values of this at intersection points will be returned.
///Straight line sl and this(RLBRep) will be projected to 2D plane of
///coordinate kind (coordinate, coordinate+1), then intersection will
///be computed. For example when sl and this are 3 dimension (x,y,z),
///and coodinate =2, 2D data (z,x) are extracted from sl and this, then
///intersection will be performed.
///**Note** MGStraight sl is treated as infinite straight line,
///even if it is finite.
MGCParam_list isect_2D(
	const MGStraight& sl,///< Straight line.
	size_t coordinate=0	///< Coordinate kind of 2D sub space.
) const;	

///Compute intersection points of 3D sub NURBS of original B-rep.
///Parameter values of thisat intersection points will be returned.
///This(RLBRep) will be projected to 3D plane of coordinate kind 
///(coordinate, coordinate+1, coordinate+2), then intersection will
///be computed. This is valid only when sdim()>=4. For example when
///pl and this are 4 dimension (x,y,z,p), and coodinate =1,
///3D data (y,z,p) are extracted from pl and this, then
///intersection will be performed.
MGCParam_list isect_3D(
	const MGPlane& pl,	///< Plane.
	size_t coordinate=0	///< Coordinate kind of 3D sub space.
) const;	

///Access to i-th element of knot
///( left-hand side version)
double& knot(size_t i){return m_line.knot(i);}			

///Access to i-th element of knot
///(right hand side version)
double knot(size_t i) const{return m_line.knot(i);}

///Returns a pointer to the knot vector data.
const double* knot_data() const{return m_line.knot_data();}		

///Returns the knot vector.
///RHS version.
const MGKnotVector& knot_vector() const{return m_line.knot_vector();}

///Returns the knot vector.
///LHS version.
MGKnotVector& knot_vector(){return m_line.knot_vector();}	

/// 自身に指定したパラメータ範囲のlimitをつける。
///Get the sub interval line of the original line.
MGRLBRep& limit(const MGInterval& itvl);

///Returns the B-coef's.
///RHS version.
const MGBPointSeq& line_bcoef() const{return m_line.line_bcoef();};

///Returns the B-coef's.
///LHS version.
MGBPointSeq& line_bcoef(){update_mark(); return m_line.line_bcoef();};

///Change direction of the line.
void negate(){m_line.negate();};

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const{return m_line.negate_param(t);};

///Return non_homogeneous B-Coefficients with weights of
///the rational B-Spline. This MGBPointSeq includes weights.
MGBPointSeq non_homogeneous_bcoef() const;

///Test if this is actually non_rational, i.e. , all of the weights are
///same values. If non_rational return true, else false.
int non_rational() const;

///Returns the order.
unsigned order() const{return m_line.order();}

/// Return ending parameter value.
double param_e() const{return m_line.param_e();};

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const{return m_line.param_normalize(t);}

/// Return starting parameter value.
double param_s() const{return m_line.param_s();};

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGRLBRep* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;

///Check if the line B-rep is planar.
///Funtion's return value is;
/// 0: Not planar, nor a point, nor straight line.
/// 1: NURBS is a point.		2: NURBS is a straight line.
/// 3: NURBS is planar.
int planar(
	MGPlane& plane	///<When Brep is not straight line nor a point,
		///< plane is returned. Even when not planar(return value is 0), plane nearest is returned.
	, MGStraight& line		///<When Brep is a line, line is returned.
	, MGPosition& point		///<When Brep is a point, point is returned.
)const;

///Change the NURBS by decreasing B-Rep dimension by ndec. This is
///an approximation of the origimal NURBS. Return value is error flag.
int reduce(
	int ndec			///<Number of B-Rep dimension to decrease 
){
	update_mark();
	return m_line.reduce(ndec);
};

///Change an original NURBS to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
MGRLBRep& refine(
		const MGKnotVector& t	///<knot vector
){
	m_line.refine(t);
	update_mark();
	return *this;
};

///ノット削除関数
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///removal knot. line_zero tolerance is used.
void remove_knot();

///Returns the space dimension.
size_t sdim() const{return m_line.sdim()-1;};

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,			///<Sweep Direction.
	double start_dist,					///<distance to start edge.
	double end_dist				///<distance to end edge.
) const;

/// 曲線のタイプをを返す。
///Return the curve type.
MGCURVE_TYPE type() const{return MGCURVE_RSPLINE;}

/// limitをはずす。
/// unlimit this line.
MGCurve& unlimit(){return *this;};

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_end(){return *this;};

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_start(){return *this;};

///Returns the knot vector.
///MGKnotVector knot_vector_real() const{return knot_vector();}

///Debug Function
public:

///IGES output function. PD126.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

std::ostream& out(std::ostream&) const;

protected:

///Compute intersection point of 1D sub NURBS of original B-rep.
///Parameter values of this at intersection points will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
)const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

std::string whoami()const{return "RLBRep";};

private:

////////////Member Data//////////////

	MGLBRep  m_line;			///< Maximum space dimension id is for weights.

///Compute the box of the whole of the curve.
///Returned is a newed object pointer.
MGBox* compute_box()const;

///Draw this by converting straight line segments.
void drawgl(
	double dl,		///<approximate line length of the straight line segments.
	double tstart, double tend	///<start and end parameter value of this.
)const;

///Compute intersection points of n-Dimensional sub NURBS and n-dimension
///plane that passes through the origin.
///Parameter values of this at intersection points will be returned.
///MGVector N is the normal vector of the plane.
MGCParam_list isect_nD(						
	const MGVector& N,		///< Normal of n-dimension plane.
	size_t dimension,		///< Number of dimension.
	size_t coordinate		///< Coordinate kind of n-dimensional sub NURBS.
) const;

///Internal function for draw_1D, _2D.
///See draw_2D for the comments of move, line, wind, and ynum.
void draw_all2D(
	int kfunc,		///<Kind of function move and line,
		///<1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], ///<window box to draw the line in.
	size_t ynum	///<Resolution of the line.
)const;

///Internal function for draw_1D, _2D.
///See draw_2D for the comments of move, line, wind, and ynum.
void draw_all1D(
	int coordinate,	///<indicates if draw_1D(>=0) or draw_2D(<0)
					///<and coordinate kind if draw_1D.
	bool t_is_x,	///<=true:t is x coordinate, and false:t is y.
	int kfunc,		///<Kind of function move and line,
		///<1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], ///<window box to draw the line in.
	size_t ynum	///<Resolution of the line.
)const;

///Split conic RLBRep at the mid point of i-th span.
///This RLBRep mmust be conic section.
MGRLBRep& split_conic(size_t i);

///2本のB表現曲線を接続する(同じ種類のとき)
///MGCurve* join(const MGCurve& crv1) const;

};

///Function to compute control point P1 and weight w1 of rational form of
///an ellipse segment. Pi and Ti are points and tangents of start and end
///for i=0,2. P is mid point of the ellipse.
///Function's output is if obtained(!=0:true) or not(=0:false).
///When obtained, =1:as finite control point, =2:as infinite.
///When T0, T2, P0, and P2 are not in one plane, function return 0.
///
///(P0,1.) (P1,w1) (P2,1.) constitute the ellipse control polygon
///of order 3 in homogeneous form.
///
///See "The NURBS Book" of W.Tiller and L.Piegl publised by Springer.
MGDECL int MGRLBRep_ellipse_weight
	(const MGPosition& P0, const MGVector& T0,
	 const MGPosition& P,
	 const MGPosition& P2, const MGVector& T2,
	 MGPosition& P1, double& w1
);

/** @} */ // end of GEO group
#endif
